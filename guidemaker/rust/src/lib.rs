use anyhow::{anyhow, Context, Result};
use bio::alphabets::dna;
use flate2::read::GzDecoder;
use hnsw_rs::prelude::*;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;
use zstd::stream::Decoder as ZstdDecoder;

/// Representation of a single parsed sequence record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqRecord {
    pub id: String,
    pub seq: Vec<u8>,
}

/// Control generation
pub mod controls;
pub use controls::{generate_flexible_negative_controls, SeedOrientation};

pub mod cfd;
pub use cfd::calculate_cfd_2bit;


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

/// Calculate precomputed left-aligned target mask
#[inline(always)]
pub fn compute_target_mask(target_len: usize) -> u64 {
    let active_bits = 2 * target_len;
    if active_bits >= 64 {
        u64::MAX
    } else {
        !((1u64 << (64 - active_bits)) - 1)
    }
}

/// Calculate exact base Hamming distance between two 2-bit left-aligned u64 sequences using masked lane-collapse
#[inline(always)]
pub fn base_hamming_distance_masked(a: u64, b: u64, mask: u64) -> u32 {
    let x = a ^ b;
    (((x | (x >> 1)) & 0x5555_5555_5555_5555) & mask).count_ones()
}

/// Get segment base count and base offset for segment s out of m segments for target_len bases
#[inline(always)]
pub fn get_segment_bounds(s: usize, m: usize, target_len: usize) -> (usize, usize) {
    let base_len = target_len / m;
    let rem = target_len % m;

    let mut start_base = 0;
    for i in 0..s {
        start_base += base_len + if i < rem { 1 } else { 0 };
    }
    let seg_bases = base_len + if s < rem { 1 } else { 0 };
    (start_base, seg_bases)
}

/// Multi-Index Hashing (MIH) CSR Table for a single segment
#[derive(Debug)]
pub struct MihSegmentTable {
    pub shift: u32,
    pub mask: u64,
    pub offsets: Vec<u32>,
    pub target_ids: Vec<u32>,
}

impl MihSegmentTable {
    pub fn build(seq_vec: &[u64], target_len: usize, s: usize, m: usize) -> Self {
        let (start_base, seg_bases) = get_segment_bounds(s, m, target_len);
        let shift = (64 - 2 * (start_base + seg_bases)) as u32;
        let mask = if 2 * seg_bases >= 64 {
            u64::MAX
        } else {
            (1u64 << (2 * seg_bases)) - 1
        };

        let num_buckets = 1usize << (2 * seg_bases);
        let mut bucket_counts = vec![0u32; num_buckets];

        for &seq in seq_vec {
            let key = ((seq >> shift) & mask) as usize;
            bucket_counts[key] += 1;
        }

        let mut offsets = vec![0u32; num_buckets + 1];
        for i in 0..num_buckets {
            offsets[i + 1] = offsets[i] + bucket_counts[i];
        }

        let mut cursor = offsets.clone();
        let mut target_ids = vec![0u32; seq_vec.len()];

        for (i, &seq) in seq_vec.iter().enumerate() {
            let key = ((seq >> shift) & mask) as usize;
            let pos = cursor[key] as usize;
            target_ids[pos] = i as u32;
            cursor[key] += 1;
        }

        MihSegmentTable {
            shift,
            mask,
            offsets,
            target_ids,
        }
    }

    #[inline(always)]
    pub fn extract_key(&self, seq: u64) -> usize {
        ((seq >> self.shift) & self.mask) as usize
    }

    #[inline(always)]
    pub fn get_bucket_targets(&self, key: usize) -> &[u32] {
        if key + 1 < self.offsets.len() {
            let start = self.offsets[key] as usize;
            let end = self.offsets[key + 1] as usize;
            &self.target_ids[start..end]
        } else {
            &[]
        }
    }

    #[inline(always)]
    pub fn get_bucket_size(&self, key: usize) -> usize {
        if key + 1 < self.offsets.len() {
            (self.offsets[key + 1] - self.offsets[key]) as usize
        } else {
            0
        }
    }
}

/// DNA Hamming Distance metric for hnsw_rs
#[derive(Clone, Copy)]
pub struct DnaHammingDistance {
    pub mask: u64,
}

impl Distance<u64> for DnaHammingDistance {
    #[inline(always)]
    fn eval(&self, va: &[u64], vb: &[u64]) -> f32 {
        base_hamming_distance_masked(va[0], vb[0], self.mask) as f32
    }
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

/// Execute Step 3 candidate reduction and top-N nearest neighbor search pipeline

// Execute Step 3 candidate reduction and top-N nearest neighbor search pipeline
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
                let score = crate::cfd::calculate_cfd_2bit(cand.seq, s, target_len);
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

// Execute Step 3 candidate reduction and top-N nearest neighbor search pipeline
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

/// Safely extract column values as a vector of String across Categorical, String, or Enum types
pub fn column_to_string_vec(df: &DataFrame, name: &str) -> Result<Vec<String>> {
    let col = df.column(name)?;
    let casted = col.cast(&DataType::String)?;
    let str_ca = casted.str()?;
    Ok(str_ca
        .into_iter()
        .map(|opt| opt.unwrap_or("").to_string())
        .collect())
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
    let mut pk_vec = Vec::with_capacity(features.len());

    for (i, f) in features.iter().enumerate() {
        pk_vec.push(i as u32);
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
        Series::new("primary_key".into(), pk_vec).into(),
        chrom_series.into(),
        Series::new("feature_start".into(), start_vec).into(),
        Series::new("feature_end".into(), end_vec).into(),
        Series::new("strand".into(), strand_vec).into(),
        Series::new("feature_id".into(), id_vec).into(),
        type_series.into(),
    ])?;

    Ok(df)
}

/// Parallel Spatial Feature-Proximity Window Filter returning matching feature primary keys Series and candidate pass mask
pub fn evaluate_spatial_filter(
    guides_df: &DataFrame,
    features_df: &DataFrame,
    before: u32,
    into: u32,
    feature_types: Option<&[String]>,
) -> Result<(Series, Vec<bool>)> {
    let num_guides = guides_df.height();

    let is_disabled = match feature_types {
        Some(ftypes) => ftypes.iter().any(|t| {
            let lower = t.trim().to_lowercase();
            lower == "disable" || lower == "none" || lower == "off"
        }),
        None => false,
    };

    if is_disabled {
        let feature_keys_series = Series::full_null("feature_keys".into(), num_guides, &DataType::List(Box::new(DataType::UInt32)));
        let init_cand_ca = guides_df.column("candidate")?.bool()?;
        let passes_vec: Vec<bool> = init_cand_ca.into_no_null_iter().collect();
        return Ok((feature_keys_series, passes_vec));
    }

    let is_all_types = match feature_types {
        None => false,
        Some(ftypes) => ftypes.iter().any(|t| t.trim().eq_ignore_ascii_case("all")),
    };

    let chrom_vec = column_to_string_vec(guides_df, "chrom")?;
    let start_ca = guides_df.column("start")?.u32()?;
    let stop_ca = guides_df.column("stop")?.u32()?;

    let feat_pks: Vec<u32> = if let Ok(col) = features_df.column("primary_key") {
        col.cast(&DataType::UInt32)?.u32()?.into_no_null_iter().collect()
    } else {
        (0..features_df.height() as u32).collect()
    };

    let feat_chrom_vec = column_to_string_vec(features_df, "chrom")?;
    let feat_start_ca = features_df.column("feature_start")?.u32()?;
    let feat_end_ca = features_df.column("feature_end")?.u32()?;
    let feat_strand_ca = features_df.column("strand")?.bool()?;
    let feat_type_vec = column_to_string_vec(features_df, "feature_type")?;

    const BIN_SHIFT: u32 = 16; // 64 kb genomic bins
    let mut chrom_bin_map: HashMap<String, HashMap<u32, Vec<(u32, u32, u32)>>> = HashMap::new();

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
        let pk = feat_pks[i];

        let (w_min, w_max) = if f_strand {
            (f_start.saturating_sub(before), f_start.saturating_add(into))
        } else {
            (f_end.saturating_sub(into), f_end.saturating_add(before))
        };

        let start_bin = w_min >> BIN_SHIFT;
        let end_bin = w_max >> BIN_SHIFT;

        let bin_map = chrom_bin_map.entry(chrom.clone()).or_default();
        for b in start_bin..=end_bin {
            bin_map.entry(b).or_default().push((w_min, w_max, pk));
        }
    }

    let chunk_size = 250_000;
    let guide_indices: Vec<usize> = (0..num_guides).collect();

    let (series_chunks, passes_chunks): (Vec<Series>, Vec<Vec<bool>>) = guide_indices
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut builder = ListPrimitiveChunkedBuilder::<UInt32Type>::new(
                "".into(),
                chunk.len(),
                chunk.len() * 2,
                DataType::UInt32,
            );
            let mut passes = Vec::with_capacity(chunk.len());

            for &i in chunk {
                let chrom = &chrom_vec[i];
                let g_start = start_ca.get(i).unwrap_or(0);
                let g_stop = stop_ca.get(i).unwrap_or(0);
                let midpoint = (g_start + g_stop) / 2;

                let mut matched_pks = Vec::new();

                if let Some(bin_map) = chrom_bin_map.get(chrom) {
                    let b = midpoint >> BIN_SHIFT;
                    if let Some(win_list) = bin_map.get(&b) {
                        for &(w_min, w_max, pk) in win_list {
                            if midpoint >= w_min && midpoint <= w_max {
                                matched_pks.push(pk);
                            }
                        }
                    }
                }

                if !matched_pks.is_empty() {
                    builder.append_slice(&matched_pks);
                    passes.push(true);
                } else {
                    builder.append_null();
                    passes.push(false);
                }
            }

            (builder.finish().into_series(), passes)
        })
        .unzip();

    let feature_keys_series = if !series_chunks.is_empty() {
        let first_series = series_chunks[0].clone();
        let mut combined_ca: ListChunked = first_series.list().unwrap().clone();
        for ch in &series_chunks[1..] {
            let list_ca: &ListChunked = ch.list().unwrap();
            combined_ca.append(list_ca).unwrap();
        }
        let mut s = combined_ca.into_series();
        s.rename("feature_keys".into());
        s
    } else {
        Series::new_empty("feature_keys".into(), &DataType::List(Box::new(DataType::UInt32)))
    };

    let passes_vec: Vec<bool> = passes_chunks.into_iter().flatten().collect();

    Ok((feature_keys_series, passes_vec))
}

/// Flat Compressed Sparse Row (CSR) Inverted Pigeonhole Index over LSR-20 region (2.55 GB RAM for 317M targets)
struct FlatCsrPigeonholeIndex {
    offsets0: Vec<u32>,
    ids0: Vec<u32>,
    offsets1: Vec<u32>,
    ids1: Vec<u32>,
}

impl FlatCsrPigeonholeIndex {
    fn build(lsr20_vec: &[u64]) -> Self {
        let num_targets = lsr20_vec.len();
        let num_buckets = 1048576; // 2^20 keys

        let mut counts0 = vec![0u32; num_buckets];
        let mut counts1 = vec![0u32; num_buckets];

        for &lsr20 in lsr20_vec {
            let k0 = ((lsr20 >> 44) as usize) & 0xFFFFF;
            let k1 = ((lsr20 >> 24) as usize) & 0xFFFFF;
            counts0[k0] += 1;
            counts1[k1] += 1;
        }

        let mut offsets0 = vec![0u32; num_buckets + 1];
        let mut offsets1 = vec![0u32; num_buckets + 1];
        for i in 0..num_buckets {
            offsets0[i + 1] = offsets0[i] + counts0[i];
            offsets1[i + 1] = offsets1[i] + counts1[i];
        }

        let mut cursor0 = offsets0.clone();
        let mut cursor1 = offsets1.clone();

        let mut ids0 = vec![0u32; num_targets];
        let mut ids1 = vec![0u32; num_targets];

        for (j, &lsr20) in lsr20_vec.iter().enumerate() {
            let k0 = ((lsr20 >> 44) as usize) & 0xFFFFF;
            let k1 = ((lsr20 >> 24) as usize) & 0xFFFFF;

            let pos0 = cursor0[k0] as usize;
            ids0[pos0] = j as u32;
            cursor0[k0] += 1;

            let pos1 = cursor1[k1] as usize;
            ids1[pos1] = j as u32;
            cursor1[k1] += 1;
        }

        FlatCsrPigeonholeIndex {
            offsets0,
            ids0,
            offsets1,
            ids1,
        }
    }

    #[inline(always)]
    fn get_postings0(&self, key0: usize) -> &[u32] {
        let start = self.offsets0[key0] as usize;
        let end = self.offsets0[key0 + 1] as usize;
        &self.ids0[start..end]
    }

    #[inline(always)]
    fn get_postings1(&self, key1: usize) -> &[u32] {
        let start = self.offsets1[key1] as usize;
        let end = self.offsets1[key1 + 1] as usize;
        &self.ids1[start..end]
    }
}

/// Multi-Threaded Parallel LSR-20 Hamming Distance <= 1 Filter across candidate guides
pub fn evaluate_lsr_hamming_filter(
    guides_df: &DataFrame,
    current_candidates: &[bool],
    target_len: usize,
    _lsr_len: usize,
    is_5prime: bool,
) -> Result<Vec<bool>> {
    let num_guides = guides_df.height();
    let seq_ca = guides_df.column("seq")?.u64()?;
    let seq_vec: Vec<u64> = seq_ca.into_no_null_iter().collect();

    let lsr_mask = compute_target_mask(20);

    let lsr20_vec: Vec<u64> = seq_vec
        .par_iter()
        .map(|&seq| {
            let key = extract_lsr_key(seq, target_len, 20, is_5prime);
            (key << 24) & lsr_mask
        })
        .collect();

    let index = FlatCsrPigeonholeIndex::build(&lsr20_vec);

    let passes: Vec<bool> = (0..num_guides)
        .into_par_iter()
        .map(|i| {
            if !current_candidates[i] {
                return false;
            }

            let cand_lsr20 = lsr20_vec[i];
            let k0 = ((cand_lsr20 >> 44) as usize) & 0xFFFFF;
            let k1 = ((cand_lsr20 >> 24) as usize) & 0xFFFFF;

            let postings0 = index.get_postings0(k0);
            for &j_u32 in postings0 {
                let j = j_u32 as usize;
                if j != i && base_hamming_distance_masked(cand_lsr20, lsr20_vec[j], lsr_mask) <= 1 {
                    return false;
                }
            }

            let postings1 = index.get_postings1(k1);
            for &j_u32 in postings1 {
                let j = j_u32 as usize;
                if j != i && base_hamming_distance_masked(cand_lsr20, lsr20_vec[j], lsr_mask) <= 1 {
                    return false;
                }
            }

            true
        })
        .collect();

    Ok(passes)
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

    let init_cand_ca = guides_df.column("candidate")?.bool()?;
    let mut current_candidates: Vec<bool> = init_cand_ca.into_no_null_iter().collect();

    // 1. Spatial interval filter FIRST
    let start_sp = Instant::now();
    let (feature_keys_series, spatial_passes) = evaluate_spatial_filter(
        guides_df,
        features_df,
        before,
        into,
        feature_types,
    )?;
    let spatial_time_sec = start_sp.elapsed().as_secs_f64();

    for i in 0..total_input_rows {
        if !spatial_passes[i] {
            current_candidates[i] = false;
        }
    }
    let rows_passing_spatial = current_candidates.iter().filter(|&&b| b).count();

    // 2. LSR-20 Hamming distance <= 1 filter SECOND
    let start_lsr = Instant::now();
    let lsr_passes = evaluate_lsr_hamming_filter(
        guides_df,
        &current_candidates,
        target_len,
        lsr_len,
        is_5prime,
    )?;
    let lsr_time_sec = start_lsr.elapsed().as_secs_f64();

    for i in 0..total_input_rows {
        if !lsr_passes[i] {
            current_candidates[i] = false;
        }
    }
    let rows_passing_lsr = current_candidates.iter().filter(|&&b| b).count();

    let final_candidate_rows = current_candidates.iter().filter(|&&b| b).count();
    let total_time_sec = start_total.elapsed().as_secs_f64();

    let mut output_df = guides_df.clone();
    let cand_series = Series::new("candidate".into(), current_candidates);
    output_df.replace("candidate", cand_series)?;
    output_df.with_column(feature_keys_series)?;

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
/// Validates that the PAM length is between 3 and 5 nt.
pub fn parse_pam_masks(pam: &str) -> Result<Vec<u8>> {
    if pam.is_empty() {
        return Err(anyhow!("PAM string cannot be empty"));
    }
    let masks: Result<Vec<u8>> = pam.chars().map(iupac_mask).collect();
    let masks = masks?;
    if masks.len() < 3 || masks.len() > 5 {
        return Err(anyhow!(
            "PAM length must be 3, 4, or 5 nt, got {}",
            masks.len()
        ));
    }
    Ok(masks)
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

/// Encode a target sequence (up to 31 nt) into a left-aligned 2-bit u64.
/// A=00, C=01, G=10, T=11.
/// Returns error if non-ATCG base encountered or length outside 1..=31.
pub fn encode_2bit_u64(seq_bytes: &[u8]) -> Result<u64> {
    if seq_bytes.is_empty() || seq_bytes.len() > 31 {
        return Err(anyhow!(
            "Target length must be between 1 and 31, got {}",
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
    let flank_len = 5usize.saturating_sub(pam_len);
    let seq_len = seq.len();
    let total_window = target_len + pam_len + flank_len;
    if seq_len < total_window {
        return hits;
    }
    let min_pos = flank_len;
    let max_pos = seq_len - pam_len - target_len;
    for pos in min_pos..=max_pos {
        if pam_matches_forward(pam_masks, seq, pos) {
            let window_start = pos - flank_len;
            let window_end = pos + pam_len + target_len;
            let window_bytes = &seq[window_start..window_end];
            if let Ok(encoded_seq) = encode_2bit_u64(window_bytes) {
                let target_start = pos + pam_len;
                let target_end = target_start + target_len;
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
    let flank_len = 5usize.saturating_sub(pam_len);
    let seq_len = seq.len();
    let total_window = target_len + pam_len + flank_len;
    if seq_len < total_window {
        return hits;
    }
    let rc_seq = dna::revcomp(seq);
    let min_rc_pos = flank_len;
    let max_rc_pos = seq_len - pam_len - target_len;
    for rc_pos in min_rc_pos..=max_rc_pos {
        if pam_matches_forward(pam_masks, &rc_seq, rc_pos) {
            let rc_window_start = rc_pos - flank_len;
            let rc_window_end = rc_pos + pam_len + target_len;
            let window_bytes = &rc_seq[rc_window_start..rc_window_end];
            if let Ok(encoded_seq) = encode_2bit_u64(window_bytes) {
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
    let flank_len = 5usize.saturating_sub(pam_len);
    let seq_len = seq.len();
    let total_window = target_len + pam_len + flank_len;
    if seq_len < total_window {
        return hits;
    }
    let min_pos = target_len;
    let max_pos = seq_len - pam_len - flank_len;
    for pos in min_pos..=max_pos {
        if pam_matches_forward(pam_masks, seq, pos) {
            let target_start = pos - target_len;
            let window_end = target_start + total_window;
            let window_bytes = &seq[target_start..window_end];
            if let Ok(encoded_seq) = encode_2bit_u64(window_bytes) {
                let start = target_start as u32;
                let stop = pos as u32;
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
    let flank_len = 5usize.saturating_sub(pam_len);
    let seq_len = seq.len();
    let total_window = target_len + pam_len + flank_len;
    if seq_len < total_window {
        return hits;
    }
    let rc_seq = dna::revcomp(seq);
    let min_rc_pos = target_len;
    let max_rc_pos = seq_len - pam_len - flank_len;
    for rc_pos in min_rc_pos..=max_rc_pos {
        if pam_matches_forward(pam_masks, &rc_seq, rc_pos) {
            let rc_target_start = rc_pos - target_len;
            let rc_window_end = rc_target_start + total_window;
            let window_bytes = &rc_seq[rc_target_start..rc_window_end];
            if let Ok(encoded_seq) = encode_2bit_u64(window_bytes) {
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
    ParquetWriter::new(file)
        .with_compression(ParquetCompression::Zstd(Some(ZstdLevel::try_new(1)?)))
        .finish(df)?;
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
    fn test_pam_length_validation() {
        assert!(parse_pam_masks("NGG").is_ok());
        assert!(parse_pam_masks("TTTN").is_ok());
        assert!(parse_pam_masks("TTTVN").is_ok());
        assert!(parse_pam_masks("NG").is_err());
        assert!(parse_pam_masks("NNGGNN").is_err());
    }

    #[test]
    fn test_encode_2bit_u64() {
        let seq = b"ACGT";
        let encoded = encode_2bit_u64(seq).unwrap();
        assert_eq!(encoded, 0x1B00_0000_0000_0000);

        assert!(encode_2bit_u64(b"ACGTN").is_err());
        let long_seq = vec![b'A'; 32];
        assert!(encode_2bit_u64(&long_seq).is_err());
    }

    #[test]
    fn test_lsr_key_extraction_5prime_and_3prime() {
        let seq = encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap();

        let key_5p = extract_lsr_key(seq, 20, 8, true);
        let expected_5p = (encode_2bit_u64(b"ACGTACGT").unwrap()) >> 48;
        assert_eq!(key_5p, expected_5p);

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
            20,
            false,
            Some(&["gene".to_string()]),
        )
        .unwrap();

        let cand_ca = filtered_df.column("candidate").unwrap().bool().unwrap();
        assert_eq!(cand_ca.get(0), Some(false));
        assert_eq!(cand_ca.get(1), Some(false));
        assert_eq!(cand_ca.get(2), Some(false));
        assert_eq!(stats.final_candidate_rows, 0);

        println!("Filtered DF columns: {:?}", filtered_df.get_column_names());
        let keys_ca = filtered_df.column("feature_keys").unwrap().list().unwrap();
        assert!(keys_ca.get_as_series(0).is_some());
        assert!(keys_ca.get_as_series(1).is_some());
        assert!(keys_ca.get_as_series(2).is_none());
    }

    #[test]
    fn test_step2_feature_types_all_and_disable() {
        let hits = vec![
            TargetHit {
                candidate: true,
                seq: encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap(),
                chrom_idx: 0,
                start: 9000,
                stop: 9020,
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
            feature_id: "exon1".to_string(),
            feature_type: "exon".to_string(),
        }];
        let features_df = build_features_dataframe(&features).unwrap();

        // 1. With "all" feature types -> passes exon feature
        let (df_all, stats_all) = execute_step2(
            &guides_df,
            &features_df,
            2000,
            500,
            20,
            20,
            false,
            Some(&["all".to_string()]),
        )
        .unwrap();
        assert_eq!(stats_all.final_candidate_rows, 1);
        assert_eq!(df_all.column("candidate").unwrap().bool().unwrap().get(0), Some(true));
        assert!(df_all.column("feature_keys").unwrap().list().unwrap().get_as_series(0).is_some());

        // 2. With "disable" feature types -> spatial filter disabled, passes based on LSR alone
        let (df_dis, stats_dis) = execute_step2(
            &guides_df,
            &features_df,
            2000,
            500,
            20,
            20,
            false,
            Some(&["disable".to_string()]),
        )
        .unwrap();
        assert_eq!(stats_dis.final_candidate_rows, 1);
        assert_eq!(df_dis.column("candidate").unwrap().bool().unwrap().get(0), Some(true));
        assert!(df_dis.column("feature_keys").unwrap().list().unwrap().get_as_series(0).is_none());
    }
}