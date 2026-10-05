use anyhow::{anyhow, Context, Result};
use polars::prelude::*;
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;
use std::time::Instant;
use crate::analysis::distance::base_hamming_distance_masked;
use crate::analysis::distance::compute_target_mask;
use crate::io::column_to_string_vec;

/// Step 3 Benchmark & Filtering statistics
#[derive(Debug, Clone)]
pub struct Step3Stats {
    pub total_candidates: usize,
    pub passed_candidates: usize,
    pub failed_candidates: usize,
    pub pass_rate: f64,
    pub index_build_time_sec: f64,
    pub search_time_sec: f64,
    pub total_time_sec: f64,
}

/// Helper container for tracking top N nearest matches during exact search
#[derive(Debug, Clone)]
pub struct TopNMatches {
    top_n: usize,
    matches: Vec<(u32, u64)>,
}

impl TopNMatches {
    fn new(top_n: usize) -> Self {
        Self {
            top_n,
            matches: Vec::with_capacity(top_n + 1),
        }
    }

    #[inline(always)]
    fn insert(&mut self, dist: u32, seq: u64) {
        if self.matches.len() < self.top_n {
            self.matches.push((dist, seq));
            self.matches.sort_unstable_by_key(|&(d, _)| d);
        } else if dist < self.matches.last().unwrap().0 {
            self.matches.pop();
            self.matches.push((dist, seq));
            self.matches.sort_unstable_by_key(|&(d, _)| d);
        }
    }

    fn into_results(self) -> (Vec<u32>, Vec<u64>) {
        let mut dists = Vec::with_capacity(self.matches.len());
        let mut seqs = Vec::with_capacity(self.matches.len());
        for (d, s) in self.matches {
            dists.push(d);
            seqs.push(s);
        }
        (dists, seqs)
    }
}

// Helper to extract 4 independent 10-bit blocks (5 bases each) from the 20 nt target spacer (bits 24..63)
#[inline(always)]
pub fn extract_pigeonhole_blocks(seq: u64) -> [usize; 4] {
    let block0 = ((seq >> 54) & 0x03FF) as usize;
    let block1 = ((seq >> 44) & 0x03FF) as usize;
    let block2 = ((seq >> 34) & 0x03FF) as usize;
    let block3 = ((seq >> 24) & 0x03FF) as usize;
    [block0, block1, block2, block3]
}

/// Representation of a candidate guide for Step 3 search
pub struct CandidateRecord {
    pub guide_df_idx: usize,
    pub seq: u64,
    pub chrom: String,
    pub start: u32,
    pub stop: u32,
    pub strand: bool,
    pub passed: AtomicBool,
    pub top_matches: Mutex<TopNMatches>,
}

/// Pigeonhole Index of candidate guides
pub struct CandidateIndex {
    pub table0: Vec<Vec<usize>>,
    pub table1: Vec<Vec<usize>>,
    pub table2: Vec<Vec<usize>>,
    pub table3: Vec<Vec<usize>>,
    pub candidates: Vec<CandidateRecord>,
}

/// Build Pigeonhole Index for candidate guides from guides_df
pub fn build_candidate_index(guides_df: &DataFrame, top_n: usize) -> Result<CandidateIndex> {
    let cand_ca = guides_df.column("candidate")?.bool()?;
    let cand_vec: Vec<bool> = cand_ca.into_no_null_iter().collect();

    let seq_ca = guides_df.column("seq")?.u64()?;
    let seq_vec: Vec<u64> = seq_ca.into_no_null_iter().collect();

    let chrom_vec = column_to_string_vec(guides_df, "chrom")?;
    let start_ca = guides_df.column("start")?.u32()?;
    let stop_ca = guides_df.column("stop")?.u32()?;
    let strand_ca = guides_df.column("strand")?.bool()?;

    let mut candidates = Vec::new();
    let mut table0 = vec![Vec::<usize>::new(); 1024];
    let mut table1 = vec![Vec::<usize>::new(); 1024];
    let mut table2 = vec![Vec::<usize>::new(); 1024];
    let mut table3 = vec![Vec::<usize>::new(); 1024];

    for (df_idx, &is_cand) in cand_vec.iter().enumerate() {
        if !is_cand {
            continue;
        }

        let cand_idx = candidates.len();
        let seq = seq_vec[df_idx];
        let chrom = chrom_vec[df_idx].clone();
        let start = start_ca.get(df_idx).unwrap_or(0);
        let stop = stop_ca.get(df_idx).unwrap_or(0);
        let strand = strand_ca.get(df_idx).unwrap_or(true);

        let blocks = extract_pigeonhole_blocks(seq);
        table0[blocks[0]].push(cand_idx);
        table1[blocks[1]].push(cand_idx);
        table2[blocks[2]].push(cand_idx);
        table3[blocks[3]].push(cand_idx);

        candidates.push(CandidateRecord {
            guide_df_idx: df_idx,
            seq,
            chrom,
            start,
            stop,
            strand,
            passed: AtomicBool::new(true),
            top_matches: Mutex::new(TopNMatches::new(top_n)),
        });
    }

    Ok(CandidateIndex {
        table0,
        table1,
        table2,
        table3,
        candidates,
    })
}

/// Scan a batch of genomic targets against the candidate guide index
pub fn scan_targets_against_candidates(
    targets_df: &DataFrame,
    index: &CandidateIndex,
    d: u32,
    target_len: usize,
) -> Result<()> {
    let num_targets = targets_df.height();
    let seq_ca = targets_df.column("seq")?.u64()?;
    let seq_vec: Vec<u64> = seq_ca.into_no_null_iter().collect();

    let chrom_vec = column_to_string_vec(targets_df, "chrom")?;
    let start_ca = targets_df.column("start")?.u32()?;
    let strand_ca = targets_df.column("strand")?.bool()?;

    let target_mask = compute_target_mask(target_len.min(20));

    (0..num_targets)
        .into_par_iter()
        .with_min_len(1024)
        .for_each(|g_idx| {
            let g_seq = seq_vec[g_idx];
            let g_chrom = &chrom_vec[g_idx];
            let g_start = start_ca.get(g_idx).unwrap_or(0);
            let g_strand = strand_ca.get(g_idx).unwrap_or(true);

            let blocks = extract_pigeonhole_blocks(g_seq);

            let process_candidate = |cand_idx: usize| {
                let cand = &index.candidates[cand_idx];

                // Skip if this genomic target is the candidate guide's own locus
                if g_strand == cand.strand && g_start == cand.start && g_chrom == &cand.chrom {
                    return;
                }

                let dist = base_hamming_distance_masked(cand.seq, g_seq, target_mask);
                if dist < d {
                    cand.passed.store(false, Ordering::Relaxed);
                } else {
                    if let Ok(mut guard) = cand.top_matches.lock() {
                        guard.insert(dist, g_seq);
                    }
                }
            };

            for &c_idx in &index.table0[blocks[0]] { process_candidate(c_idx); }
            for &c_idx in &index.table1[blocks[1]] { process_candidate(c_idx); }
            for &c_idx in &index.table2[blocks[2]] { process_candidate(c_idx); }
            for &c_idx in &index.table3[blocks[3]] { process_candidate(c_idx); }
        });

    Ok(())
}

/// Compile final Polars DataFrame for passing candidate guides from CandidateIndex
pub fn build_candidate_results_dataframe(
    guides_df: &DataFrame,
    index: CandidateIndex,
    top_n: usize,
    target_len: usize,
) -> Result<DataFrame> {
    let mut passed_df_indices = Vec::new();
    let mut passed_dists_nested = Vec::new();
    let mut passed_seqs_nested = Vec::new();
    let mut passed_cfds_nested = Vec::new();

    for cand in index.candidates {
        if cand.passed.load(Ordering::Relaxed) {
            passed_df_indices.push(cand.guide_df_idx);

            let top_matches = cand.top_matches.into_inner().unwrap();
            let (dists, seqs) = top_matches.into_results();

            let mut cfds: Vec<f32> = Vec::with_capacity(seqs.len());
            for &s in &seqs {
                let score = crate::analysis::cfd::calculate_cfd_2bit(cand.seq, s, target_len);
                cfds.push(score);
            }

            let mut padded_dists: Vec<u32> = dists;
            let mut padded_seqs: Vec<u64> = seqs;
            let mut padded_cfds: Vec<f32> = cfds;
            while padded_dists.len() < top_n {
                padded_dists.push(0);
                padded_seqs.push(0);
                padded_cfds.push(1.0);
            }

            passed_dists_nested.push(padded_dists);
            passed_seqs_nested.push(padded_seqs);
            passed_cfds_nested.push(padded_cfds);
        }
    }

    let mut passed_mask = vec![false; guides_df.height()];
    for &df_idx in &passed_df_indices {
        passed_mask[df_idx] = true;
    }

    let mask_ca = ChunkedArray::<BooleanType>::from_slice("cand_mask".into(), &passed_mask);
    let filtered_df = guides_df.filter(&mask_ca)?;

    let dist_series_vec: Vec<Series> = passed_dists_nested
        .into_iter()
        .map(|v| Series::new("".into(), v))
        .collect();
    let nn_dist_series = Series::new("nn_dist".into(), dist_series_vec);

    let seq_series_vec: Vec<Series> = passed_seqs_nested
        .into_iter()
        .map(|v| Series::new("".into(), v))
        .collect();
    let nn_seq_series = Series::new("nn_seq".into(), seq_series_vec);

    let cfd_series_vec: Vec<Series> = passed_cfds_nested
        .into_iter()
        .map(|v| Series::new("".into(), v))
        .collect();
    let nn_cfd_series = Series::new("nn_cfd".into(), cfd_series_vec);

    let mut final_df = filtered_df;
    final_df.with_column(nn_dist_series)?;
    final_df.with_column(nn_seq_series)?;
    final_df.with_column(nn_cfd_series)?;

    Ok(final_df)
}

/// Execute Step 3 candidate reduction and top-N nearest neighbor search pipeline
pub fn execute_step3(
    guides_df: &DataFrame,
    d: u32,
    top_n: usize,
    _method: &str,
    target_len: usize,
) -> Result<(DataFrame, Step3Stats)> {
    let start_total = Instant::now();

    let start_idx = Instant::now();
    let index = build_candidate_index(guides_df, top_n)?;
    let index_build_time_sec = start_idx.elapsed().as_secs_f64();

    let start_search = Instant::now();
    scan_targets_against_candidates(guides_df, &index, d, target_len)?;
    let search_time_sec = start_search.elapsed().as_secs_f64();

    let total_candidates = index.candidates.len();
    let final_df = build_candidate_results_dataframe(guides_df, index, top_n, target_len)?;

    let passed_candidates = final_df.height();
    let failed_candidates = total_candidates.saturating_sub(passed_candidates);
    let pass_rate = if total_candidates > 0 {
        (passed_candidates as f64 / total_candidates as f64) * 100.0
    } else {
        0.0
    };

    let total_time_sec = start_total.elapsed().as_secs_f64();

    let stats = Step3Stats {
        total_candidates,
        passed_candidates,
        failed_candidates,
        pass_rate,
        index_build_time_sec,
        search_time_sec,
        total_time_sec,
    };

    Ok((final_df, stats))
}
