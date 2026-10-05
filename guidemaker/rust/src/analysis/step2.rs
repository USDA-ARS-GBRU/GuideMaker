use anyhow::Result;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Instant;
use crate::analysis::distance::base_hamming_distance_masked;
use crate::io::column_to_string_vec;

/// Step 2 Benchmark statistics
#[derive(Debug, Clone)]
pub struct Step2Stats {
    pub total_input_rows: usize,
    pub rows_passing_spatial: usize,
    pub rows_passing_lsr: usize,
    pub final_candidate_rows: usize,
    pub spatial_time_sec: f64,
    pub lsr_time_sec: f64,
    pub total_time_sec: f64,
}

/// Extract LSR key from 2-bit u64 sequence depending on orientation
#[inline(always)]
pub fn extract_lsr_key(seq: u64, target_len: usize, lsr_len: usize, is_5prime: bool) -> u64 {
    if is_5prime {
        seq >> (64 - 2 * lsr_len)
    } else {
        let mask = if lsr_len >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * lsr_len)) - 1
        };
        (seq >> (64 - 2 * target_len)) & mask
    }
}

/// Parallel Spatial Feature-Proximity Window Filter evaluated ONLY on active candidate rows
pub fn evaluate_spatial_filter(
    guides_df: &DataFrame,
    features_df: &DataFrame,
    current_candidates: &[bool],
    before: u32,
    into: u32,
    feature_types: Option<&[String]>,
) -> Result<Vec<bool>> {
    let num_guides = guides_df.height();

    if let Some(ftypes) = feature_types {
        if ftypes.iter().any(|t| {
            let lower = t.trim().to_lowercase();
            lower == "disable" || lower == "none" || lower == "off"
        }) {
            return Ok(current_candidates.to_vec());
        }
    }

    let is_all_types = match feature_types {
        None => false,
        Some(ftypes) => ftypes.iter().any(|t| t.trim().eq_ignore_ascii_case("all")),
    };

    let chrom_vec = column_to_string_vec(guides_df, "chrom")?;
    let start_ca = guides_df.column("start")?.u32()?;
    let stop_ca = guides_df.column("stop")?.u32()?;

    let feat_chrom_vec = column_to_string_vec(features_df, "chrom")?;
    let feat_start_ca = features_df.column("feature_start")?.u32()?;
    let feat_end_ca = features_df.column("feature_end")?.u32()?;
    let feat_strand_ca = features_df.column("strand")?.bool()?;
    let feat_type_vec = column_to_string_vec(features_df, "feature_type")?;

    let mut chrom_tss_map: HashMap<String, Vec<u32>> = HashMap::new();

    for i in 0..features_df.height() {
        if !is_all_types {
            if let Some(ftypes) = feature_types {
                let ftype_str = &feat_type_vec[i];
                if !ftypes.iter().any(|t| t.trim().eq_ignore_ascii_case(ftype_str)) {
                    continue;
                }
            }
        }

        let chrom = &feat_chrom_vec[i];
        let f_start = feat_start_ca.get(i).unwrap_or(0);
        let f_end = feat_end_ca.get(i).unwrap_or(0);
        let f_strand = feat_strand_ca.get(i).unwrap_or(true);

        let tss = if f_strand { f_start } else { f_end };
        chrom_tss_map.entry(chrom.clone()).or_default().push(tss);
    }

    for tss_vec in chrom_tss_map.values_mut() {
        tss_vec.sort_unstable();
    }

    let passes: Vec<bool> = (0..num_guides)
        .into_par_iter()
        .map(|i| {
            if !current_candidates[i] {
                return false;
            }

            let chrom = &chrom_vec[i];
            let g_start = start_ca.get(i).unwrap_or(0);
            let g_stop = stop_ca.get(i).unwrap_or(0);
            let midpoint = (g_start + g_stop) / 2;

            if let Some(tss_list) = chrom_tss_map.get(chrom) {
                let upper = midpoint.saturating_add(before);

                let idx = if midpoint < into {
                    0
                } else {
                    let lower = midpoint - into;
                    tss_list.partition_point(|&x| x <= lower)
                };

                if idx < tss_list.len() && tss_list[idx] < upper {
                    return true;
                }
            }

            false
        })
        .collect();

    Ok(passes)
}

/// Multi-Threaded Parallel LSR Uniqueness Filter across active candidate guides
pub fn evaluate_lsr_uniqueness(
    seq_slice: &[u64],
    current_candidates: &[bool],
    target_len: usize,
    lsr_len: usize,
    is_5prime: bool,
) -> Vec<bool> {
    if lsr_len <= 14 {
        let table_size = 1usize << (2 * lsr_len);
        let counts: Vec<AtomicU32> = (0..table_size).map(|_| AtomicU32::new(0)).collect();

        seq_slice
            .par_iter()
            .enumerate()
            .for_each(|(i, &seq)| {
                if current_candidates[i] {
                    let key = extract_lsr_key(seq, target_len, lsr_len, is_5prime) as usize;
                    counts[key].fetch_add(1, Ordering::Relaxed);
                }
            });

        seq_slice
            .par_iter()
            .enumerate()
            .map(|(i, &seq)| {
                if current_candidates[i] {
                    let key = extract_lsr_key(seq, target_len, lsr_len, is_5prime) as usize;
                    counts[key].load(Ordering::Relaxed) == 1
                } else {
                    false
                }
            })
            .collect()
    } else {
        let counts: HashMap<u64, u32> = seq_slice
            .par_iter()
            .enumerate()
            .fold(
                || HashMap::new(),
                |mut acc, (i, &seq)| {
                    if current_candidates[i] {
                        let key = extract_lsr_key(seq, target_len, lsr_len, is_5prime);
                        *acc.entry(key).or_insert(0) += 1;
                    }
                    acc
                },
            )
            .reduce(
                || HashMap::new(),
                |mut map1, map2| {
                    for (k, v) in map2 {
                        *map1.entry(k).or_insert(0) += v;
                    }
                    map1
                },
            );

        seq_slice
            .par_iter()
            .enumerate()
            .map(|(i, &seq)| {
                if current_candidates[i] {
                    let key = extract_lsr_key(seq, target_len, lsr_len, is_5prime);
                    counts.get(&key) == Some(&1)
                } else {
                    false
                }
            })
            .collect()
    }
}

/// Execute Step-2 Pipeline and collect benchmark statistics
pub fn execute_step2(
    guides_df: &DataFrame,
    features_df: &DataFrame,
    before: u32,
    into: u32,
    target_len: usize,
    lsr_len: usize,
    is_5prime: bool,
    feature_types: Option<&[String]>,
) -> Result<(DataFrame, Step2Stats)> {
    let start_total = Instant::now();
    let total_input_rows = guides_df.height();

    let seq_ca = guides_df.column("seq")?.u64()?;
    let seq_vec: Vec<u64> = seq_ca.into_no_null_iter().collect();

    let mut current_candidates = vec![true; total_input_rows];

    // Evaluate LSR uniqueness first
    let start_lsr = Instant::now();
    let lsr_passes = evaluate_lsr_uniqueness(&seq_vec, &current_candidates, target_len, lsr_len, is_5prime);
    let lsr_time_sec = start_lsr.elapsed().as_secs_f64();

    for i in 0..total_input_rows {
        if !lsr_passes[i] {
            current_candidates[i] = false;
        }
    }
    let rows_passing_lsr = current_candidates.iter().filter(|&&b| b).count();

    // Evaluate Spatial filter strictly on LSR candidate rows
    let start_sp = Instant::now();
    let spatial_passes = evaluate_spatial_filter(guides_df, features_df, &current_candidates, before, into, feature_types)?;
    let spatial_time_sec = start_sp.elapsed().as_secs_f64();

    for i in 0..total_input_rows {
        if !spatial_passes[i] {
            current_candidates[i] = false;
        }
    }
    let rows_passing_spatial = current_candidates.iter().filter(|&&b| b).count();

    let final_candidate_rows = current_candidates.iter().filter(|&&b| b).count();
    let total_time_sec = start_total.elapsed().as_secs_f64();

    let mut output_df = guides_df.clone();
    let cand_series = Series::new("candidate".into(), current_candidates);
    output_df.replace("candidate", cand_series)?;

    let stats = Step2Stats {
        total_input_rows,
        rows_passing_spatial,
        rows_passing_lsr,
        final_candidate_rows,
        spatial_time_sec,
        lsr_time_sec,
        total_time_sec,
    };

    Ok((output_df, stats))
}
