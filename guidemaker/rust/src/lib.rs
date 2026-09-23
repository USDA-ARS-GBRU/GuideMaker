use anyhow::{anyhow, Context, Result};
use bio::alphabets::dna;
use flate2::read::GzDecoder;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Instant;
use zstd::stream::Decoder as ZstdDecoder;

/// Representation of a single parsed sequence record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqRecord {
    pub id: String,
    pub seq: Vec<u8>,
}

/// Representation of a single candidate hit
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TargetHit {
    pub candidate: bool,
    pub seq: u64,
    pub chrom_idx: u32,
    pub start: u32,
    pub stop: u32,
    pub strand: bool,
}

/// Representation of a genomic feature record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureRecord {
    pub chrom: String,
    pub feature_start: u32,
    pub feature_end: u32,
    pub strand: bool,
    pub feature_id: String,
    pub feature_type: String,
}

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

/// Transparently open a file with optional gzip or zstd decompression
pub fn open_compressed_reader(path: &Path) -> Result<Box<dyn Read + Send>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open sequence file at {:?}", path))?;
    let mut buf_reader = BufReader::new(file);

    let mut magic = [0u8; 4];
    let n = buf_reader.read(&mut magic)?;

    let combined_reader = std::io::Cursor::new(magic[..n].to_vec()).chain(buf_reader);

    if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        let gz_decoder = GzDecoder::new(combined_reader);
        Ok(Box::new(gz_decoder))
    } else if n >= 4 && magic[0] == 0x28 && magic[1] == 0xb5 && magic[2] == 0x2f && magic[3] == 0xfd {
        let zstd_decoder = ZstdDecoder::new(combined_reader)
            .context("Failed to initialize Zstd decoder")?;
        Ok(Box::new(zstd_decoder))
    } else {
        Ok(Box::new(combined_reader))
    }
}

/// Parse FASTA records from reader
pub fn parse_fasta_records<R: Read>(reader: R) -> Result<Vec<SeqRecord>> {
    let fasta_reader = bio::io::fasta::Reader::new(reader);
    let mut records = Vec::new();
    for rec_res in fasta_reader.records() {
        let rec = rec_res.context("Failed to parse FASTA record")?;
        records.push(SeqRecord {
            id: rec.id().to_string(),
            seq: rec.seq().to_ascii_uppercase(),
        });
    }
    Ok(records)
}

/// Parse GenBank records from reader
pub fn parse_genbank_records<R: Read>(reader: R) -> Result<Vec<SeqRecord>> {
    let buf_reader = BufReader::new(reader);
    let mut records = Vec::new();

    let mut current_id = String::new();
    let mut current_seq = Vec::new();
    let mut in_origin = false;

    for line_res in buf_reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();

        if trimmed.starts_with("LOCUS") {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() > 1 {
                current_id = parts[1].to_string();
            } else {
                current_id = "unknown_contig".to_string();
            }
            current_seq.clear();
            in_origin = false;
        } else if trimmed.starts_with("ACCESSION") && current_id.is_empty() {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() > 1 {
                current_id = parts[1].to_string();
            }
        } else if trimmed.starts_with("VERSION") && current_id.is_empty() {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() > 1 {
                current_id = parts[1].to_string();
            }
        } else if trimmed.starts_with("ORIGIN") {
            in_origin = true;
        } else if trimmed == "//" {
            if !current_id.is_empty() && !current_seq.is_empty() {
                records.push(SeqRecord {
                    id: std::mem::take(&mut current_id),
                    seq: std::mem::take(&mut current_seq),
                });
            } else {
                current_id.clear();
                current_seq.clear();
            }
            in_origin = false;
        } else if in_origin {
            for b in line.bytes() {
                if b.is_ascii_alphabetic() {
                    current_seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }

    if !current_id.is_empty() && !current_seq.is_empty() {
        records.push(SeqRecord {
            id: current_id,
            seq: current_seq,
        });
    }

    Ok(records)
}

/// Read FASTA or GenBank sequence records transparently
pub fn read_sequence_records(path: &Path) -> Result<Vec<SeqRecord>> {
    let mut reader = open_compressed_reader(path)?;

    let mut header_buf = [0u8; 1024];
    let n = reader.read(&mut header_buf)?;

    let peek_str = String::from_utf8_lossy(&header_buf[..n]);

    let combined_reader = std::io::Cursor::new(header_buf[..n].to_vec()).chain(reader);

    if peek_str.trim_start().starts_with('>') {
        parse_fasta_records(combined_reader)
    } else if peek_str.contains("LOCUS") || peek_str.contains("ORIGIN") || peek_str.contains("FEATURES") {
        parse_genbank_records(combined_reader)
    } else {
        parse_fasta_records(combined_reader)
    }
}

/// Parse GFF/GTF features from reader
pub fn parse_gff_gtf_features<R: Read>(reader: R) -> Result<Vec<FeatureRecord>> {
    let buf_reader = BufReader::new(reader);
    let mut features = Vec::new();

    for line_res in buf_reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }

        let chrom = cols[0].to_string();
        let feature_type = cols[2].to_string();
        let start_1based: u32 = match cols[3].parse() {
            Ok(val) => val,
            Err(_) => continue,
        };
        let end_1based: u32 = match cols[4].parse() {
            Ok(val) => val,
            Err(_) => continue,
        };

        let feature_start = start_1based.saturating_sub(1);
        let feature_end = end_1based;
        let strand = cols[6] != "-";

        let attrs = cols[8];
        let feature_id = extract_gff_gtf_id(attrs, &chrom, feature_start, feature_end, &feature_type);

        features.push(FeatureRecord {
            chrom,
            feature_start,
            feature_end,
            strand,
            feature_id,
            feature_type,
        });
    }

    Ok(features)
}

fn extract_gff_gtf_id(
    attrs: &str,
    chrom: &str,
    start: u32,
    end: u32,
    ftype: &str,
) -> String {
    let mut found_id = String::new();

    for part in attrs.split(';') {
        let trimmed = part.trim();
        if let Some((key, val)) = trimmed.split_once('=') {
            let k = key.trim();
            let clean_val = val.trim_matches('"').trim();
            if !clean_val.is_empty() {
                if k.eq_ignore_ascii_case("locus_tag") || k.eq_ignore_ascii_case("ID") {
                    return clean_val.to_string();
                } else if (k.eq_ignore_ascii_case("gene_id") || k.eq_ignore_ascii_case("Name") || k.eq_ignore_ascii_case("transcript_id"))
                    && found_id.is_empty()
                {
                    found_id = clean_val.to_string();
                }
            }
        }
    }

    if !found_id.is_empty() {
        return found_id;
    }

    for part in attrs.split(';') {
        let trimmed = part.trim();
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.len() >= 2 {
            let k = tokens[0];
            let clean_val = tokens[1].trim_matches('"').trim();
            if !clean_val.is_empty() {
                if k.eq_ignore_ascii_case("locus_tag") || k.eq_ignore_ascii_case("gene_id") {
                    return clean_val.to_string();
                } else if k.eq_ignore_ascii_case("transcript_id") && found_id.is_empty() {
                    found_id = clean_val.to_string();
                }
            }
        }
    }

    if !found_id.is_empty() {
        return found_id;
    }

    format!("{}_{}_{}_{}", chrom, start, end, ftype)
}

/// Parse GenBank FEATURES section into FeatureRecord list
pub fn parse_genbank_features<R: Read>(reader: R) -> Result<Vec<FeatureRecord>> {
    let buf_reader = BufReader::new(reader);
    let mut features = Vec::new();

    let mut current_chrom = "unknown_contig".to_string();
    let mut in_features = false;
    let mut current_type = String::new();
    let mut current_loc = String::new();
    let mut current_id = String::new();

    for line_res in buf_reader.lines() {
        let line = line_res?;

        if line.starts_with("LOCUS") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() > 1 {
                current_chrom = parts[1].to_string();
            }
            in_features = false;
        } else if line.starts_with("FEATURES") {
            in_features = true;
        } else if line.starts_with("ORIGIN") || line.starts_with("//") {
            if !current_type.is_empty() {
                if let Some((f_start, f_end, strand)) = parse_genbank_location(&current_loc) {
                    let fid = if !current_id.is_empty() {
                        std::mem::take(&mut current_id)
                    } else {
                        format!("{}_{}_{}_{}", current_chrom, f_start, f_end, current_type)
                    };
                    features.push(FeatureRecord {
                        chrom: current_chrom.clone(),
                        feature_start: f_start,
                        feature_end: f_end,
                        strand,
                        feature_id: fid,
                        feature_type: std::mem::take(&mut current_type),
                    });
                } else {
                    current_type.clear();
                    current_id.clear();
                }
            }
            in_features = false;
        } else if in_features {
            if line.len() >= 21 && !line[5..21].trim().is_empty() && !line[5..21].contains('/') {
                if !current_type.is_empty() {
                    if let Some((f_start, f_end, strand)) = parse_genbank_location(&current_loc) {
                        let fid = if !current_id.is_empty() {
                            std::mem::take(&mut current_id)
                        } else {
                            format!("{}_{}_{}_{}", current_chrom, f_start, f_end, current_type)
                        };
                        features.push(FeatureRecord {
                            chrom: current_chrom.clone(),
                            feature_start: f_start,
                            feature_end: f_end,
                            strand,
                            feature_id: fid,
                            feature_type: std::mem::take(&mut current_type),
                        });
                    }
                    current_type.clear();
                    current_id.clear();
                    current_loc.clear();
                }

                current_type = line[5..21].trim().to_string();
                current_loc = line[21..].trim().to_string();
            } else if !current_type.is_empty() {
                let trimmed = line.trim();
                if trimmed.starts_with('/') {
                    if let Some((key, val)) = trimmed[1..].split_once('=') {
                        let k = key.trim();
                        let clean_v = val.trim().trim_matches('"').trim();
                        let is_locus_tag = k.eq_ignore_ascii_case("locus_tag");
                        let is_gene_id = k.eq_ignore_ascii_case("gene")
                            || k.eq_ignore_ascii_case("protein_id")
                            || k.eq_ignore_ascii_case("db_xref")
                            || k.eq_ignore_ascii_case("ID");

                        if is_locus_tag || (is_gene_id && current_id.is_empty()) {
                            current_id = clean_v.to_string();
                        }
                    }
                } else if !current_loc.is_empty() && !trimmed.starts_with('/') {
                    current_loc.push_str(trimmed);
                }
            }
        }
    }

    Ok(features)
}

fn parse_genbank_location(loc_str: &str) -> Option<(u32, u32, bool)> {
    let strand = !loc_str.contains("complement");

    let mut numbers = Vec::new();
    let mut current_num = String::new();

    for b in loc_str.bytes() {
        if b.is_ascii_digit() {
            current_num.push(b as char);
        } else if !current_num.is_empty() {
            if let Ok(n) = current_num.parse::<u32>() {
                numbers.push(n);
            }
            current_num.clear();
        }
    }
    if !current_num.is_empty() {
        if let Ok(n) = current_num.parse::<u32>() {
            numbers.push(n);
        }
    }

    if numbers.is_empty() {
        return None;
    }

    let min_pos = *numbers.iter().min()?;
    let max_pos = *numbers.iter().max()?;

    if min_pos == 0 {
        return None;
    }

    let feature_start = min_pos - 1;
    let feature_end = max_pos;

    Some((feature_start, feature_end, strand))
}

/// Read FeatureRecord list transparently from GFF/GTF or GenBank file
pub fn read_feature_records(path: &Path) -> Result<Vec<FeatureRecord>> {
    let mut reader = open_compressed_reader(path)?;

    let mut header_buf = [0u8; 1024];
    let n = reader.read(&mut header_buf)?;

    let peek_str = String::from_utf8_lossy(&header_buf[..n]);

    let combined_reader = std::io::Cursor::new(header_buf[..n].to_vec()).chain(reader);

    if peek_str.contains("LOCUS") || peek_str.contains("FEATURES") || peek_str.contains("ORIGIN") {
        parse_genbank_features(combined_reader)
    } else {
        parse_gff_gtf_features(combined_reader)
    }
}

/// Build Polars features DataFrame from FeatureRecord list
pub fn build_features_dataframe(features: &[FeatureRecord]) -> Result<DataFrame> {
    let mut chrom_vec = Vec::with_capacity(features.len());
    let mut start_vec = Vec::with_capacity(features.len());
    let mut end_vec = Vec::with_capacity(features.len());
    let mut strand_vec = Vec::with_capacity(features.len());
    let mut id_vec = Vec::with_capacity(features.len());
    let mut type_vec = Vec::with_capacity(features.len());

    for f in features {
        chrom_vec.push(f.chrom.as_str());
        start_vec.push(f.feature_start);
        end_vec.push(f.feature_end);
        strand_vec.push(f.strand);
        id_vec.push(f.feature_id.as_str());
        type_vec.push(f.feature_type.as_str());
    }

    let chrom_series = Series::new("chrom".into(), chrom_vec)
        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;
    let type_series = Series::new("feature_type".into(), type_vec)
        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;

    let df = DataFrame::new(vec![
        chrom_series.into(),
        Series::new("feature_start".into(), start_vec).into(),
        Series::new("feature_end".into(), end_vec).into(),
        Series::new("strand".into(), strand_vec).into(),
        Series::new("feature_id".into(), id_vec).into(),
        type_series.into(),
    ])?;

    Ok(df)
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
    let chrom_str_col = guides_df.column("chrom")?.cast(&DataType::String)?;
    let chrom_ca = chrom_str_col.str()?;
    let start_ca = guides_df.column("start")?.u32()?;
    let stop_ca = guides_df.column("stop")?.u32()?;

    let feat_chrom_str_col = features_df.column("chrom")?.cast(&DataType::String)?;
    let feat_chrom_ca = feat_chrom_str_col.str()?;
    let feat_start_ca = features_df.column("feature_start")?.u32()?;
    let feat_end_ca = features_df.column("feature_end")?.u32()?;
    let feat_strand_ca = features_df.column("strand")?.bool()?;

    let feat_type_str_col = features_df.column("feature_type")?.cast(&DataType::String)?;
    let feat_type_ca = feat_type_str_col.str()?;

    let mut chrom_tss_map: HashMap<String, Vec<u32>> = HashMap::new();

    for i in 0..features_df.height() {
        if let Some(ftypes) = feature_types {
            let ftype_str = feat_type_ca.get(i).unwrap_or("");
            if !ftypes.iter().any(|t| t.eq_ignore_ascii_case(ftype_str)) {
                continue;
            }
        }

        let chrom = feat_chrom_ca.get(i).unwrap_or("");
        let f_start = feat_start_ca.get(i).unwrap_or(0);
        let f_end = feat_end_ca.get(i).unwrap_or(0);
        let f_strand = feat_strand_ca.get(i).unwrap_or(true);

        let tss = if f_strand { f_start } else { f_end };
        chrom_tss_map.entry(chrom.to_string()).or_default().push(tss);
    }

    for tss_vec in chrom_tss_map.values_mut() {
        tss_vec.sort_unstable();
    }

    let num_guides = guides_df.height();

    // Rayon parallel evaluation across candidate guides
    let passes: Vec<bool> = (0..num_guides)
        .into_par_iter()
        .map(|i| {
            if !current_candidates[i] {
                return false;
            }

            let chrom = chrom_ca.get(i).unwrap_or("");
            let g_start = start_ca.get(i).unwrap_or(0);
            let g_stop = stop_ca.get(i).unwrap_or(0);
            let midpoint = (g_start + g_stop) / 2;

            if let Some(tss_list) = chrom_tss_map.get(chrom) {
                let lower = midpoint.saturating_sub(into);
                let upper = midpoint.saturating_add(before);

                let idx = tss_list.partition_point(|&x| x <= lower);
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

        // Rayon parallel count pass
        seq_slice
            .par_iter()
            .enumerate()
            .for_each(|(i, &seq)| {
                if current_candidates[i] {
                    let key = extract_lsr_key(seq, target_len, lsr_len, is_5prime) as usize;
                    counts[key].fetch_add(1, Ordering::Relaxed);
                }
            });

        // Rayon parallel uniqueness evaluation pass
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
        // Fallback parallel map-reduce for large lsr_len
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
    fast_filter_first: bool,
) -> Result<(DataFrame, Step2Stats)> {
    let start_total = Instant::now();
    let total_input_rows = guides_df.height();

    let seq_ca = guides_df.column("seq")?.u64()?;
    let seq_vec: Vec<u64> = seq_ca.into_no_null_iter().collect();

    let mut current_candidates = vec![true; total_input_rows];

    let rows_passing_spatial;
    let rows_passing_lsr;
    let spatial_time_sec;
    let lsr_time_sec;

    if fast_filter_first {
        // Strategy A: Spatial filter first, then LSR uniqueness
        let start_sp = Instant::now();
        let spatial_passes = evaluate_spatial_filter(guides_df, features_df, &current_candidates, before, into, feature_types)?;
        spatial_time_sec = start_sp.elapsed().as_secs_f64();

        for i in 0..total_input_rows {
            if !spatial_passes[i] {
                current_candidates[i] = false;
            }
        }
        rows_passing_spatial = current_candidates.iter().filter(|&&b| b).count();

        let start_lsr = Instant::now();
        let lsr_passes = evaluate_lsr_uniqueness(&seq_vec, &current_candidates, target_len, lsr_len, is_5prime);
        lsr_time_sec = start_lsr.elapsed().as_secs_f64();

        for i in 0..total_input_rows {
            if !lsr_passes[i] {
                current_candidates[i] = false;
            }
        }
        rows_passing_lsr = current_candidates.iter().filter(|&&b| b).count();
    } else {
        // Strategy B: LSR uniqueness first, then Spatial filter (evaluated ONLY on LSR candidates)
        let start_lsr = Instant::now();
        let lsr_passes = evaluate_lsr_uniqueness(&seq_vec, &current_candidates, target_len, lsr_len, is_5prime);
        lsr_time_sec = start_lsr.elapsed().as_secs_f64();

        for i in 0..total_input_rows {
            if !lsr_passes[i] {
                current_candidates[i] = false;
            }
        }
        rows_passing_lsr = current_candidates.iter().filter(|&&b| b).count();

        let start_sp = Instant::now();
        let spatial_passes = evaluate_spatial_filter(guides_df, features_df, &current_candidates, before, into, feature_types)?;
        spatial_time_sec = start_sp.elapsed().as_secs_f64();

        for i in 0..total_input_rows {
            if !spatial_passes[i] {
                current_candidates[i] = false;
            }
        }
        rows_passing_spatial = current_candidates.iter().filter(|&&b| b).count();
    }

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

/// Convert an IUPAC character to its bitmask.
/// A=1, C=2, G=4, T=8.
pub fn iupac_mask(c: char) -> Result<u8> {
    match c.to_ascii_uppercase() {
        'A' => Ok(1),
        'C' => Ok(2),
        'G' => Ok(4),
        'T' => Ok(8),
        'M' => Ok(1 | 2),
        'R' => Ok(1 | 4),
        'W' => Ok(1 | 8),
        'S' => Ok(2 | 4),
        'Y' => Ok(2 | 8),
        'K' => Ok(4 | 8),
        'V' => Ok(1 | 2 | 4),
        'H' => Ok(1 | 2 | 8),
        'D' => Ok(1 | 4 | 8),
        'B' => Ok(2 | 4 | 8),
        'X' | 'N' => Ok(1 | 2 | 4 | 8),
        other => Err(anyhow!("Invalid IUPAC character: '{}'", other)),
    }
}

/// Convert a PAM string into a vector of per-base IUPAC bitmasks.
pub fn parse_pam_masks(pam: &str) -> Result<Vec<u8>> {
    if pam.is_empty() {
        return Err(anyhow!("PAM string cannot be empty"));
    }
    pam.chars().map(iupac_mask).collect()
}

/// Check if PAM matches the sequence at position `pos`.
pub fn pam_matches_forward(pam_masks: &[u8], seq_bytes: &[u8], pos: usize) -> bool {
    if pos + pam_masks.len() > seq_bytes.len() {
        return false;
    }
    for (i, &pam_mask) in pam_masks.iter().enumerate() {
        let seq_char = seq_bytes[pos + i] as char;
        let seq_mask = match iupac_mask(seq_char) {
            Ok(m) => m,
            Err(_) => return false,
        };
        if (pam_mask & seq_mask) == 0 {
            return false;
        }
    }
    true
}

/// Encode a target sequence (up to 26 nt) into a left-aligned 2-bit u64.
/// A=00, C=01, G=10, T=11.
/// Returns error if non-ATCG base encountered or length outside 1..=26.
pub fn encode_2bit_u64(seq_bytes: &[u8]) -> Result<u64> {
    if seq_bytes.is_empty() || seq_bytes.len() > 26 {
        return Err(anyhow!(
            "Target length must be between 1 and 26, got {}",
            seq_bytes.len()
        ));
    }
    let mut val: u64 = 0;
    for (i, &b) in seq_bytes.iter().enumerate() {
        let code: u64 = match b.to_ascii_uppercase() {
            b'A' => 0b00,
            b'C' => 0b01,
            b'G' => 0b10,
            b'T' => 0b11,
            other => return Err(anyhow!("Non-ATCG base encountered in target: '{}'", other as char)),
        };
        let shift = 64 - 2 * (i + 1);
        val |= code << shift;
    }
    Ok(val)
}

/// Reverse complement wrapper using rust-bio
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    dna::revcomp(seq)
}

/// Search 5prime orientation on forward strand: 5'-[PAM][TARGET]-3'
pub fn search_5prime_forward(
    chrom_idx: u32,
    seq: &[u8],
    pam_masks: &[u8],
    target_len: usize,
) -> Vec<TargetHit> {
    let mut hits = Vec::new();
    let pam_len = pam_masks.len();
    let seq_len = seq.len();
    if seq_len < pam_len + target_len {
        return hits;
    }
    let max_pos = seq_len - pam_len - target_len;
    for pos in 0..=max_pos {
        if pam_matches_forward(pam_masks, seq, pos) {
            let target_start = pos + pam_len;
            let target_end = target_start + target_len;
            let target_bytes = &seq[target_start..target_end];
            if let Ok(encoded_seq) = encode_2bit_u64(target_bytes) {
                let start = target_start as u32;
                let stop = target_end as u32;
                if start < stop && (stop as usize) <= seq_len {
                    hits.push(TargetHit {
                        candidate: true,
                        seq: encoded_seq,
                        chrom_idx,
                        start,
                        stop,
                        strand: true,
                    });
                }
            }
        }
    }
    hits
}

/// Search 5prime orientation on reverse strand
pub fn search_5prime_reverse(
    chrom_idx: u32,
    seq: &[u8],
    pam_masks: &[u8],
    target_len: usize,
) -> Vec<TargetHit> {
    let mut hits = Vec::new();
    let pam_len = pam_masks.len();
    let seq_len = seq.len();
    if seq_len < pam_len + target_len {
        return hits;
    }
    let rc_seq = dna::revcomp(seq);
    let max_rc_pos = seq_len - pam_len - target_len;
    for rc_pos in 0..=max_rc_pos {
        if pam_matches_forward(pam_masks, &rc_seq, rc_pos) {
            let rc_target_start = rc_pos + pam_len;
            let rc_target_end = rc_target_start + target_len;
            let target_bytes = &rc_seq[rc_target_start..rc_target_end];
            if let Ok(encoded_seq) = encode_2bit_u64(target_bytes) {
                let start_fwd = seq_len - (rc_pos + pam_len + target_len);
                let stop_fwd = seq_len - (rc_pos + pam_len);
                if start_fwd < stop_fwd && stop_fwd <= seq_len {
                    hits.push(TargetHit {
                        candidate: true,
                        seq: encoded_seq,
                        chrom_idx,
                        start: start_fwd as u32,
                        stop: stop_fwd as u32,
                        strand: false,
                    });
                }
            }
        }
    }
    hits
}

/// Search 3prime orientation on forward strand: 5'-[TARGET][PAM]-3'
pub fn search_3prime_forward(
    chrom_idx: u32,
    seq: &[u8],
    pam_masks: &[u8],
    target_len: usize,
) -> Vec<TargetHit> {
    let mut hits = Vec::new();
    let pam_len = pam_masks.len();
    let seq_len = seq.len();
    if seq_len < pam_len + target_len {
        return hits;
    }
    let min_pos = target_len;
    let max_pos = seq_len - pam_len;
    for pos in min_pos..=max_pos {
        if pam_matches_forward(pam_masks, seq, pos) {
            let target_start = pos - target_len;
            let target_end = pos;
            let target_bytes = &seq[target_start..target_end];
            if let Ok(encoded_seq) = encode_2bit_u64(target_bytes) {
                let start = target_start as u32;
                let stop = target_end as u32;
                if start < stop && (stop as usize) <= seq_len {
                    hits.push(TargetHit {
                        candidate: true,
                        seq: encoded_seq,
                        chrom_idx,
                        start,
                        stop,
                        strand: true,
                    });
                }
            }
        }
    }
    hits
}

/// Search 3prime orientation on reverse strand
pub fn search_3prime_reverse(
    chrom_idx: u32,
    seq: &[u8],
    pam_masks: &[u8],
    target_len: usize,
) -> Vec<TargetHit> {
    let mut hits = Vec::new();
    let pam_len = pam_masks.len();
    let seq_len = seq.len();
    if seq_len < pam_len + target_len {
        return hits;
    }
    let rc_seq = dna::revcomp(seq);
    let min_rc_pos = target_len;
    let max_rc_pos = seq_len - pam_len;
    for rc_pos in min_rc_pos..=max_rc_pos {
        if pam_matches_forward(pam_masks, &rc_seq, rc_pos) {
            let rc_target_start = rc_pos - target_len;
            let rc_target_end = rc_pos;
            let target_bytes = &rc_seq[rc_target_start..rc_target_end];
            if let Ok(encoded_seq) = encode_2bit_u64(target_bytes) {
                let start_fwd = seq_len - rc_pos;
                let stop_fwd = seq_len - (rc_pos - target_len);
                if start_fwd < stop_fwd && stop_fwd <= seq_len {
                    hits.push(TargetHit {
                        candidate: true,
                        seq: encoded_seq,
                        chrom_idx,
                        start: start_fwd as u32,
                        stop: stop_fwd as u32,
                        strand: false,
                    });
                }
            }
        }
    }
    hits
}

/// Build Polars DataFrame from hits
pub fn build_dataframe(hits: &[TargetHit], chrom_names: &[String]) -> Result<DataFrame> {
    let mut candidate_vec = Vec::with_capacity(hits.len());
    let mut seq_vec = Vec::with_capacity(hits.len());
    let mut chrom_vec = Vec::with_capacity(hits.len());
    let mut start_vec = Vec::with_capacity(hits.len());
    let mut stop_vec = Vec::with_capacity(hits.len());
    let mut strand_vec = Vec::with_capacity(hits.len());

    for hit in hits {
        candidate_vec.push(hit.candidate);
        seq_vec.push(hit.seq);
        let name = chrom_names.get(hit.chrom_idx as usize).map(|s| s.as_str()).unwrap_or("unknown");
        chrom_vec.push(name);
        start_vec.push(hit.start);
        stop_vec.push(hit.stop);
        strand_vec.push(hit.strand);
    }

    let chrom_series = Series::new("chrom".into(), chrom_vec)
        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;

    let df = DataFrame::new(vec![
        Series::new("candidate".into(), candidate_vec).into(),
        Series::new("seq".into(), seq_vec).into(),
        chrom_series.into(),
        Series::new("start".into(), start_vec).into(),
        Series::new("stop".into(), stop_vec).into(),
        Series::new("strand".into(), strand_vec).into(),
    ])?;

    Ok(df)
}

/// Write Polars DataFrame to CSV
pub fn write_csv(df: &mut DataFrame, path: &Path) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create CSV file at {:?}", path))?;
    CsvWriter::new(file).finish(df)?;
    Ok(())
}

/// Write Polars DataFrame to Parquet
pub fn write_parquet(df: &mut DataFrame, path: &Path) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create Parquet file at {:?}", path))?;
    ParquetWriter::new(file).finish(df)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iupac_mask() {
        assert_eq!(iupac_mask('A').unwrap(), 1);
        assert_eq!(iupac_mask('C').unwrap(), 2);
        assert_eq!(iupac_mask('G').unwrap(), 4);
        assert_eq!(iupac_mask('T').unwrap(), 8);
        assert_eq!(iupac_mask('N').unwrap(), 15);
        assert_eq!(iupac_mask('M').unwrap(), 3);
        assert!(iupac_mask('Z').is_err());
    }

    #[test]
    fn test_encode_2bit_u64() {
        let seq = b"ACGT";
        let encoded = encode_2bit_u64(seq).unwrap();
        assert_eq!(encoded, 0x1B00_0000_0000_0000);

        assert!(encode_2bit_u64(b"ACGTN").is_err());
        let long_seq = vec![b'A'; 27];
        assert!(encode_2bit_u64(&long_seq).is_err());
    }

    #[test]
    fn test_lsr_key_extraction_5prime_and_3prime() {
        // 20-mer sequence: ACGTACGTACGTACGTACGT
        let seq = encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap();

        // 5prime orientation: 5' (left) end 8 nt = "ACGTACGT"
        let key_5p = extract_lsr_key(seq, 20, 8, true);
        let expected_5p = (encode_2bit_u64(b"ACGTACGT").unwrap()) >> 48;
        assert_eq!(key_5p, expected_5p);

        // 3prime orientation: 3' (right) end 8 nt = "ACGTACGT"
        let key_3p = extract_lsr_key(seq, 20, 8, false);
        assert_eq!(key_3p, expected_5p);
    }

    #[test]
    fn test_step2_spatial_and_lsr_filtering() {
        let hits = vec![
            TargetHit {
                candidate: true,
                seq: encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap(),
                chrom_idx: 0,
                start: 9000,
                stop: 9020,
                strand: true,
            },
            TargetHit {
                candidate: true,
                seq: encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap(),
                chrom_idx: 0,
                start: 9500,
                stop: 9520,
                strand: true,
            },
            TargetHit {
                candidate: true,
                seq: encode_2bit_u64(b"TGCATGCATGCATGCATGCA").unwrap(),
                chrom_idx: 0,
                start: 50000,
                stop: 50020,
                strand: true,
            },
        ];
        let chrom_names = vec!["chr1".to_string()];
        let guides_df = build_dataframe(&hits, &chrom_names).unwrap();

        let features = vec![FeatureRecord {
            chrom: "chr1".to_string(),
            feature_start: 10000,
            feature_end: 15000,
            strand: true,
            feature_id: "gene1".to_string(),
            feature_type: "gene".to_string(),
        }];
        let features_df = build_features_dataframe(&features).unwrap();

        let (filtered_df, stats) = execute_step2(
            &guides_df,
            &features_df,
            2000,
            500,
            20,
            8,
            false, // 3prime orientation
            None,
            false, // LSR first (default)
        )
        .unwrap();

        let cand_ca = filtered_df.column("candidate").unwrap().bool().unwrap();
        assert_eq!(cand_ca.get(0), Some(false));
        assert_eq!(cand_ca.get(1), Some(false));
        assert_eq!(cand_ca.get(2), Some(false));
        assert_eq!(stats.final_candidate_rows, 0);
    }
}
