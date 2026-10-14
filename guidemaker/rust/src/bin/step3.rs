use anyhow::{Context, Result};
use clap::Parser;
use fxhash::FxHashMap;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;
use sysinfo::System;

// Arrow / Parquet
use arrow_array::{
    Array, BooleanArray, DictionaryArray, Float32Array, LargeStringArray, RecordBatch, StringArray,
    UInt32Array, UInt64Array,
};
use arrow_array::types::UInt32Type;
use arrow_schema::{DataType, Field, Schema};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

// CFD (YAML-loaded; mismatch-only evaluator)
use guidemaker_scan::cfd::{CfdMatrix, load_cfd_matrix, cfd_evaluate_lowbits};

/// Fixed LSR length: 20 nt (closest to PAM)
const LSR_LEN: usize = 20;

/// MIH radius r = 4 => m = r + 1 partitions
const MIH_MISMATCH_RADIUS: u32 = 4;
const MIH_PARTITIONS: usize = (MIH_MISMATCH_RADIUS as usize) + 1; // 5
const MIH_SLICE_LEN_NT: usize = 4; // 4 nt per slice
const MIH_SLICE_OFFSETS: [usize; MIH_PARTITIONS] = [0, 4, 8, 12, 16];

/// All 10 slice pairs (a < b)
const SLICE_PAIRS: [(usize, usize); 10] = [
    (0, 1), (0, 2), (0, 3), (0, 4),
    (1, 2), (1, 3), (1, 4),
    (2, 3), (2, 4),
    (3, 4),
];

/// Fallback union budget (max IDs to scan when intersection finds no hit)
const FALLBACK_BUDGET: usize = 10_000;

/// Key cardinality per slice (4-nt → 8-bit key)
const KEY_CARDINALITY: usize = 256;

/// CLI args
#[derive(Parser, Debug)]
#[command(author, version, about = "Step-3 (Arrow) Exact MIH & Max-CFD (3′) with YAML-loaded CFD + Parallel Indexing + 3-slice bitset refine")]
pub struct Step3Args {
    /// Path to input Parquet file (unified target/guide dataset)
    #[arg(long, visible_alias = "guides", visible_alias = "input")]
    pub input: PathBuf,

    /// LSR length for Hamming & CFD evaluation (default 20 nt)
    #[arg(long, default_value_t = 20)]
    pub lsr_len: usize,

    /// Arrow batch size (rows per RecordBatch)
    #[arg(long, default_value_t = 1_048_576)]
    pub batch_size: usize,

    /// Whether to keep non-candidate rows in output (we output candidates only)
    #[arg(long, default_value_t = false)]
    pub keep_noncandidates: bool,

    /// Max candidate depth per slice bucket (skip buckets exceeding this)
    #[arg(long, default_value_t = usize::MAX)]
    pub max_bucket_depth: usize,

    /// Output Parquet file path
    #[arg(long)]
    pub out: PathBuf,

    /// Number of worker threads (default: all available CPUs)
    #[arg(long)]
    pub threads: Option<usize>,

    /// Path to a YAML file containing one or more CFD matrices
    #[arg(long, default_value = "cfd_matrices.yaml")]
    pub cfd_yaml: PathBuf,

    /// Name of the CFD matrix to use (e.g., "SpCas9_v1")
    #[arg(long, default_value = "SpCas9_v1")]
    pub cfd_name: String,
}

/// Candidate entry stored in the MIH index
#[derive(Clone)]
struct Candidate {
    seq2b_low: u64, // 20 nt packed into low bits (2-bit per nt) — 3′ case assumed
    chrom: String,
    start: u32,
    strand: bool,
    // best match fields (updated during search)
    best_score: f32,
    best_match_seq2b_low: u64,
    best_match_chrom: String,
    best_match_start: u32,
    best_match_strand: bool,
}

/// Per-bucket structure: sorted IDs + bitset membership
struct BucketIndex {
    ids: Vec<usize>,  // sorted candidate IDs
    bits: Vec<u64>,   // bitset [(N+63)/64] words
}

/// 5 slice buckets: key -> BucketIndex
type SliceBucketMap = FxHashMap<u32, BucketIndex>;

/// Precomputed 2-slice pair table: [key_a][key_b] -> Vec<u32> (sorted intersection)
struct PairTable {
    inter: Vec<Vec<u32>>, // length = KEY_CARDINALITY * KEY_CARDINALITY, row-major
}
impl PairTable {
    fn new() -> Self {
        Self {
            inter: (0..(KEY_CARDINALITY * KEY_CARDINALITY))
                .map(|_| Vec::<u32>::new())
                .collect(),
        }
    }
    #[inline]
    fn row_mut(&mut self, key_a: usize) -> &mut [Vec<u32>] {
        let start = key_a * KEY_CARDINALITY;
        &mut self.inter[start..start + KEY_CARDINALITY]
    }
    #[inline]
    fn get(&self, key_a: u32, key_b: u32) -> &Vec<u32> {
        &self.inter[(key_a as usize) * KEY_CARDINALITY + (key_b as usize)]
    }
}

/// MIH index: 5 slice maps, 10 pair tables, and candidates
struct MIHIndex {
    buckets: [SliceBucketMap; MIH_PARTITIONS],
    pair_tables: [PairTable; SLICE_PAIRS.len()],
    candidates: Vec<Candidate>, // index by cand_id
    total_candidates: usize,
    // Diagnostics (global)
    total_bucket_probes: u64,
    total_candidates_deduped: u64,
    total_cfd_evaluations: u64,
}
impl MIHIndex {
    fn new() -> Self {
        Self {
            buckets: [
                FxHashMap::default(),
                FxHashMap::default(),
                FxHashMap::default(),
                FxHashMap::default(),
                FxHashMap::default(),
            ],
            pair_tables: [
                PairTable::new(), PairTable::new(), PairTable::new(), PairTable::new(), PairTable::new(),
                PairTable::new(), PairTable::new(), PairTable::new(), PairTable::new(), PairTable::new(),
            ],
            candidates: Vec::new(),
            total_candidates: 0,
            total_bucket_probes: 0,
            total_candidates_deduped: 0,
            total_cfd_evaluations: 0,
        }
    }
}

/// === Low-bit ops (3′ case only) ===

#[inline]
fn repack_high_to_lowbits(x: u64) -> u64 {
    let mut y = 0u64;
    for pos in 0..LSR_LEN {
        let src_shift = 64 - 2 * (pos + 1);
        let base = (x >> src_shift) & 0b11;
        y |= base << (2 * pos);
    }
    y
}

#[inline]
fn slice_key_4nt_lowbits(seq2b_low: u64, off: usize) -> u32 {
    let shift = 2 * off;
    ((seq2b_low >> shift) & 0xFF) as u32
}

#[inline]
fn mismatch_mask_20_lowbits(tgt_low: u64, cand_low: u64) -> u32 {
    let mut x = tgt_low ^ cand_low;
    x |= x >> 1;
    ((x & 0x5555_5555_5555_5555u64) & ((1u64 << 20) - 1)) as u32
}

/// === Arrow helpers ===

fn col_idx(batch: &RecordBatch, name: &str) -> usize {
    batch
        .schema()
        .fields()
        .iter()
        .position(|f| f.name() == name)
        .unwrap_or_else(|| panic!("column '{}' not found", name))
}
#[inline]
fn col_u64(batch: &RecordBatch, idx: usize) -> &UInt64Array {
    batch.column(idx).as_any().downcast_ref::<UInt64Array>().expect("u64 col")
}
#[inline]
fn col_u32(batch: &RecordBatch, idx: usize) -> &UInt32Array {
    batch.column(idx).as_any().downcast_ref::<UInt32Array>().expect("u32 col")
}
#[inline]
fn col_bool(batch: &RecordBatch, idx: usize) -> &BooleanArray {
    batch.column(idx).as_any().downcast_ref::<BooleanArray>().expect("bool col")
}
#[inline]
fn get_str_at(batch: &RecordBatch, idx: usize, i: usize) -> &str {
    let col = batch.column(idx);
    match col.data_type() {
        DataType::Utf8 => {
            let arr = col.as_any().downcast_ref::<StringArray>().expect("Utf8");
            if arr.is_valid(i) { arr.value(i) } else { "" }
        }
        DataType::LargeUtf8 => {
            let arr = col.as_any().downcast_ref::<LargeStringArray>().expect("LargeUtf8");
            if arr.is_valid(i) { arr.value(i) } else { "" }
        }
        DataType::Dictionary(_, v) if **v == DataType::Utf8 => {
            let dict = col.as_any().downcast_ref::<DictionaryArray<UInt32Type>>().expect("Dict<UInt32,Utf8>");
            let keys = dict.keys();
            if keys.is_null(i) { "" } else {
                let key = keys.value(i) as usize;
                let values = dict.values().as_any().downcast_ref::<StringArray>().expect("dict Utf8 vals");
                values.value(key)
            }
        }
        _ => "",
    }
}

/// === Index construction (Arrow streaming) ===
/// (unchanged logic; only CFD scoring later is swapped to the new evaluator)

fn build_mih_index_arrow(input: &PathBuf, batch_size: usize, max_bucket_depth: usize) -> Result<MIHIndex> {
    let file = File::open(input).with_context(|| format!("Failed to open {:?}", input))?;
    let mut reader = ParquetRecordBatchReaderBuilder::try_new(file)?
        .with_batch_size(batch_size)
        .build()?;

    let mut index = MIHIndex::new();

    // Pass 1: stream candidates & fill bucket IDs
    while let Some(batch_res) = reader.next() {
        let batch = batch_res?;

        let i_candidate = col_idx(&batch, "candidate");
        let i_seq = col_idx(&batch, "seq");
        let i_chrom = col_idx(&batch, "chrom");
        let i_start = col_idx(&batch, "start");
        let i_strand = col_idx(&batch, "strand");

        let cand = col_bool(&batch, i_candidate);
        let seq = col_u64(&batch, i_seq);
        let start = col_u32(&batch, i_start);
        let strand = col_bool(&batch, i_strand);

        for r in 0..batch.num_rows() {
            if cand.is_valid(r) && cand.value(r) {
                let seq_low = repack_high_to_lowbits(seq.value(r));
                let chrom_s = get_str_at(&batch, i_chrom, r).to_string();
                let start_u32 = start.value(r);
                let strand_b = strand.value(r);

                let cid = index.candidates.len();
                index.candidates.push(Candidate {
                    seq2b_low: seq_low,
                    chrom: chrom_s,
                    start: start_u32,
                    strand: strand_b,
                    best_score: 0.0,
                    best_match_seq2b_low: 0,
                    best_match_chrom: String::new(),
                    best_match_start: 0,
                    best_match_strand: false,
                });

                for (bi, &off) in MIH_SLICE_OFFSETS.iter().enumerate() {
                    let key = slice_key_4nt_lowbits(seq_low, off);
                    index
                        .buckets[bi]
                        .entry(key)
                        .or_insert_with(|| BucketIndex { ids: Vec::new(), bits: Vec::new() })
                        .ids
                        .push(cid); // sorted by construction
                }
            }
        }
    }

    index.total_candidates = index.candidates.len();

    // Optional pruning
    if max_bucket_depth < usize::MAX {
        for bi in 0..MIH_PARTITIONS {
            index.buckets[bi].retain(|_, b| b.ids.len() <= max_bucket_depth);
        }
    }

    // Pass 2: precompute 2-slice intersections (parallel)
    index
        .pair_tables
        .par_iter_mut()
        .zip(SLICE_PAIRS.par_iter().copied())
        .for_each(|(table, (a, b))| {
            table
                .inter
                .par_chunks_mut(KEY_CARDINALITY) // one row per key_a
                .enumerate()
                .for_each(|(key_a, row)| {
                    let key_a_u32 = key_a as u32;
                    let ids_a = index.buckets[a]
                        .get(&key_a_u32)
                        .map(|bb| bb.ids.as_slice())
                        .unwrap_or(&[]);
                    if ids_a.is_empty() {
                        for out in row.iter_mut() { out.clear(); }
                        return;
                    }
                    for key_b in 0..KEY_CARDINALITY {
                        let key_b_u32 = key_b as u32;
                        let ids_b = index.buckets[b]
                            .get(&key_b_u32)
                            .map(|bb| bb.ids.as_slice())
                            .unwrap_or(&[]);
                        let out = &mut row[key_b];
                        out.clear();
                        if ids_b.is_empty() { continue; }
                        // linear merge intersection (both sorted)
                        let (mut i, mut j) = (0usize, 0usize);
                        while i < ids_a.len() && j < ids_b.len() {
                            let va = ids_a[i];
                            let vb = ids_b[j];
                            if va == vb {
                                out.push(va as u32);
                                i += 1; j += 1;
                            } else if va < vb {
                                i += 1;
                            } else {
                                j += 1;
                            }
                        }
                    }
                });
        });

    // Pass 3: build bitsets (fast 3-slice refine)
    let words = (index.total_candidates + 63) / 64;
    for bi in 0..MIH_PARTITIONS {
        index.buckets[bi].par_iter_mut().for_each(|(_k, bucket)| {
            bucket.bits = vec![0u64; words];
            for &id in bucket.ids.iter() {
                let w = id >> 6;
                let b = id & 63;
                unsafe {
                    *bucket.bits.get_unchecked_mut(w) |= 1u64 << b;
                }
            }
        });
    }

    Ok(index)
}

/// Budgeted union (sorted, dedup)
#[inline]
fn union_sorted_with_budget(lists: &[&[usize]], budget: usize) -> Vec<usize> {
    let mut out: Vec<usize> = Vec::with_capacity(budget.min(4096));
    for &lst in lists {
        let mut tmp: Vec<usize> = Vec::with_capacity(out.len() + lst.len());
        let (mut i, mut j) = (0usize, 0usize);
        let mut last: Option<usize> = None;
        while (i < out.len() || j < lst.len()) && tmp.len() < budget {
            let pick_out = if j >= lst.len() { true } else if i >= out.len() { false } else { out[i] <= lst[j] };
            let val = if pick_out { let v = out[i]; i += 1; v } else { let v = lst[j]; j += 1; v };
            if last.map_or(true, |prev| prev != val) {
                tmp.push(val); last = Some(val);
            }
        }
        out = tmp;
        if out.len() >= budget { break; }
    }
    out
}

/// Batch stats
#[derive(Default, Clone)]
struct BatchStats {
    rows: u64,
    sum_bucket_depths: [u64; MIH_PARTITIONS],
    sum_intersection2_size: u64, // |I12|
    sum_intersection3_size: u64, // |I123|
    sum_fallback_scanned: u64,
    sum_unique_scanned: u64,
    sum_hamming_pass: u64,
    sum_cfd_evals: u64,
    sum_collect_ns: u128,
    sum_hamming_ns: u128,
}

/// Choose best pair (smallest combined bucket sizes)
#[inline]
fn choose_best_pair(
    buckets: &[SliceBucketMap; MIH_PARTITIONS],
    keys: [u32; MIH_PARTITIONS],
) -> (usize, usize, usize) {
    let mut best_idx = 0usize;
    let mut best_pair = SLICE_PAIRS[0];
    let mut best_score = usize::MAX;
    for (idx, &(a, b)) in SLICE_PAIRS.iter().enumerate() {
        let la = buckets[a].get(&keys[a]).map(|bb| bb.ids.len()).unwrap_or(0);
        let lb = buckets[b].get(&keys[b]).map(|bb| bb.ids.len()).unwrap_or(0);
        let score = la + lb;
        if score < best_score {
            best_score = score;
            best_idx = idx;
            best_pair = (a, b);
        }
    }
    (best_idx, best_pair.0, best_pair.1)
}

/// Choose smallest third slice among remaining three
#[inline]
fn choose_third_slice(
    buckets: &[SliceBucketMap; MIH_PARTITIONS],
    keys: [u32; MIH_PARTITIONS],
    a: usize,
    b: usize,
) -> usize {
    let mut best_c = usize::MAX;
    let mut best_len = usize::MAX;
    for c in 0..MIH_PARTITIONS {
        if c == a || c == b { continue; }
        let len = buckets[c].get(&keys[c]).map(|bb| bb.ids.len()).unwrap_or(0);
        if len < best_len { best_len = len; best_c = c; }
    }
    best_c
}

/// Bitset refine: I12 ∩ third slice → I123
#[inline]
fn refine_i12_with_third_bitset(
    i12: &[u32],
    third_bits: &[u64],
    out: &mut Vec<usize>,
) {
    out.clear();
    out.reserve(i12.len());
    for &id_u32 in i12.iter() {
        let id = id_u32 as usize;
        if (third_bits[id >> 6] & (1u64 << (id & 63))) != 0 {
            out.push(id);
        }
    }
}

/// Process one RecordBatch of targets (uses YAML CFD matrix)
fn process_targets_batch_against_index(
    batch: &RecordBatch,
    index: &mut MIHIndex,
    lsr_len: usize,
    radius: u32,
    cfd_matrix: &Arc<CfdMatrix>,
) -> Result<BatchStats> {
    let idx_seq = col_idx(batch, "seq");
    let idx_chrom = col_idx(batch, "chrom");
    let idx_start = col_idx(batch, "start");
    let idx_strand = col_idx(batch, "strand");

    let seq = col_u64(batch, idx_seq);
    let start = col_u32(batch, idx_start);
    let strand = col_bool(batch, idx_strand);

    let n = batch.num_rows();

    let stats = (0..n)
        .into_par_iter()
        .map(|row| {
            let seq_hi = seq.value(row);
            let chrom_s = get_str_at(batch, idx_chrom, row).to_string();
            let chrom_arc: Arc<str> = Arc::from(chrom_s);
            let start_u32 = start.value(row);
            let strand_b = strand.value(row);

            // Low-bit target
            let tgt_seq_low = repack_high_to_lowbits(seq_hi);

            // Candidate collection
            let t_collect_start = Instant::now();

            let mut keys: [u32; MIH_PARTITIONS] = [0; MIH_PARTITIONS];
            for (bi, &off) in MIH_SLICE_OFFSETS.iter().enumerate() {
                keys[bi] = slice_key_4nt_lowbits(tgt_seq_low, off);
            }

            let mut bucket_depths: [u64; MIH_PARTITIONS] = [0; MIH_PARTITIONS];
            for bi in 0..MIH_PARTITIONS {
                bucket_depths[bi] = index.buckets[bi].get(&keys[bi]).map(|bb| bb.ids.len() as u64).unwrap_or(0);
            }

            let (pair_idx, a, b) = choose_best_pair(&index.buckets, keys);
            let c = choose_third_slice(&index.buckets, keys, a, b);

            // Precomputed I12
            let i12_u32 = index.pair_tables[pair_idx].get(keys[a], keys[b]);
            let i12_size = i12_u32.len();

            // I123 via bitset
            let mut i123: Vec<usize> = Vec::with_capacity(8);
            if i12_size > 0 {
                if let Some(bc) = index.buckets[c].get(&keys[c]) {
                    refine_i12_with_third_bitset(i12_u32, bc.bits.as_slice(), &mut i123);
                }
            }

            let mut fallback_scanned = 0usize;
            let collect_ns;

            // Hamming + CFD
            let mut row_updates: Vec<(usize, f32, u64, Arc<str>, u32, bool)> =
                Vec::with_capacity(i123.len().max(i12_size).max(32));
            let mut hamming_pass: u64 = 0;
            let mut cfd_eval: u64 = 0;

            if !i123.is_empty() {
                collect_ns = t_collect_start.elapsed().as_nanos();

                let t_hamming_start = Instant::now();
                for &cid in i123.iter() {
                    let cseq = index.candidates[cid].seq2b_low;
                    let mism = mismatch_mask_20_lowbits(tgt_seq_low, cseq);
                    let dist = mism.count_ones();
                    if dist <= radius {
                        hamming_pass += 1;
                        let score = cfd_evaluate_lowbits(cfd_matrix, tgt_seq_low, cseq);
                        cfd_eval += 1;
                        row_updates.push((cid, score, tgt_seq_low, chrom_arc.clone(), start_u32, strand_b));
                    }
                }
                let hamming_ns = t_hamming_start.elapsed().as_nanos();

                (
                    row_updates, bucket_depths, i12_size as u64, i123.len() as u64,
                    fallback_scanned as u64, i123.len() as u64, hamming_pass, cfd_eval, collect_ns, hamming_ns
                )
            } else if i12_size > 0 {
                collect_ns = t_collect_start.elapsed().as_nanos();

                let t_hamming_start = Instant::now();
                for &id_u32 in i12_u32.iter() {
                    let cid = id_u32 as usize;
                    let cseq = index.candidates[cid].seq2b_low;
                    let mism = mismatch_mask_20_lowbits(tgt_seq_low, cseq);
                    let dist = mism.count_ones();
                    if dist <= radius {
                        hamming_pass += 1;
                        let score = cfd_evaluate_lowbits(cfd_matrix, tgt_seq_low, cseq);
                        cfd_eval += 1;
                        row_updates.push((cid, score, tgt_seq_low, chrom_arc.clone(), start_u32, strand_b));
                    }
                }
                let hamming_ns = t_hamming_start.elapsed().as_nanos();

                (
                    row_updates, bucket_depths, i12_size as u64, 0u64,
                    fallback_scanned as u64, i12_size as u64, hamming_pass, cfd_eval, collect_ns, hamming_ns
                )
            } else {
                // Fallback union
                let t_collect_mid = Instant::now();
                let mut blists: Vec<(&[usize], usize)> = Vec::with_capacity(5);
                for bi in 0..MIH_PARTITIONS {
                    let lst = index.buckets[bi].get(&keys[bi]).map(|bb| bb.ids.as_slice()).unwrap_or(&[]);
                    blists.push((lst, lst.len()));
                }
                blists.sort_by_key(|(_, len)| *len);
                let lists: Vec<&[usize]> = blists.iter().map(|(lst, _)| *lst).collect();
                let to_scan = union_sorted_with_budget(&lists, FALLBACK_BUDGET);
                fallback_scanned = to_scan.len();
                collect_ns = t_collect_mid.elapsed().as_nanos();

                let t_hamming_start = Instant::now();
                for &cid in to_scan.iter() {
                    let cseq = index.candidates[cid].seq2b_low;
                    let mism = mismatch_mask_20_lowbits(tgt_seq_low, cseq);
                    let dist = mism.count_ones();
                    if dist <= radius {
                        hamming_pass += 1;
                        let score = cfd_evaluate_lowbits(cfd_matrix, tgt_seq_low, cseq);
                        cfd_eval += 1;
                        row_updates.push((cid, score, tgt_seq_low, chrom_arc.clone(), start_u32, strand_b));
                    }
                }
                let hamming_ns = t_hamming_start.elapsed().as_nanos();

                (
                    row_updates, bucket_depths, 0u64, 0u64,
                    fallback_scanned as u64, fallback_scanned as u64, hamming_pass, cfd_eval, collect_ns, hamming_ns
                )
            }
        })
        .reduce(
            || {
                (
                    Vec::<(usize, f32, u64, Arc<str>, u32, bool)>::new(),
                    [0u64; MIH_PARTITIONS],
                    0u64, 0u64, 0u64, 0u64, 0u64, 0u64, 0u128, 0u128
                )
            },
            |mut a, b| {
                let (
                    bvec, b_bucket_depths, b_i12, b_i123, b_fallback, b_unique_scanned, b_hpass, b_cfd, b_collect_ns, b_hamming_ns
                ) = b;
                a.0.extend(bvec.into_iter());
                for bi in 0..MIH_PARTITIONS { a.1[bi] += b_bucket_depths[bi]; }
                a.2 += b_i12;
                a.3 += b_i123;
                a.4 += b_fallback;
                a.5 += b_unique_scanned;
                a.6 += b_hpass;
                a.7 += b_cfd;
                a.8 += b_collect_ns;
                a.9 += b_hamming_ns;
                a
            },
        );

    // Apply best-score updates
    let mut batch_stats = BatchStats::default();
    batch_stats.rows = n as u64;
    batch_stats.sum_bucket_depths = stats.1;
    batch_stats.sum_intersection2_size = stats.2;
    batch_stats.sum_intersection3_size = stats.3;
    batch_stats.sum_fallback_scanned = stats.4;
    batch_stats.sum_unique_scanned = stats.5;
    batch_stats.sum_hamming_pass = stats.6;
    batch_stats.sum_cfd_evals = stats.7;
    batch_stats.sum_collect_ns = stats.8;
    batch_stats.sum_hamming_ns = stats.9;

    for (cid, score, seq_low, chrom_arc, start, strand) in stats.0.into_iter() {
        let c = &mut index.candidates[cid];
        if score > c.best_score {
            c.best_score = score;
            c.best_match_seq2b_low = seq_low;
            c.best_match_chrom = chrom_arc.to_string();
            c.best_match_start = start;
            c.best_match_strand = strand;
        }
    }

    index.total_cfd_evaluations += batch_stats.sum_cfd_evals;
    index.total_candidates_deduped += batch_stats.sum_unique_scanned;
    index.total_bucket_probes += (n as u64) * (MIH_PARTITIONS as u64);

    Ok(batch_stats)
}

/// Write results (candidates + best_* metrics)
fn write_results_parquet(out_path: &PathBuf, index: &MIHIndex) -> Result<()> {
    let seq_arr = UInt64Array::from(
        index.candidates.iter().map(|c| c.seq2b_low).collect::<Vec<u64>>()
    );
    let chrom_arr = StringArray::from(
        index.candidates.iter().map(|c| c.best_match_chrom.as_str()).collect::<Vec<&str>>()
    );
    let start_arr = UInt32Array::from(
        index.candidates.iter().map(|c| c.best_match_start).collect::<Vec<u32>>()
    );
    let strand_arr = BooleanArray::from(
        index.candidates.iter().map(|c| c.best_match_strand).collect::<Vec<bool>>()
    );

    let best_cfd_arr = Float32Array::from(
        index.candidates.iter().map(|c| c.best_score).collect::<Vec<f32>>()
    );
    let best_match_seq_arr = UInt64Array::from(
        index.candidates.iter().map(|c| c.best_match_seq2b_low).collect::<Vec<u64>>()
    );

    let schema = Arc::new(Schema::new(vec![
        Field::new("seq2b_low", DataType::UInt64, false),
        Field::new("best_match_chrom", DataType::Utf8, false),
        Field::new("best_match_start", DataType::UInt32, false),
        Field::new("best_match_strand", DataType::Boolean, false),
        Field::new("best_cfd", DataType::Float32, false),
        Field::new("best_match_seq2b_low", DataType::UInt64, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(seq_arr),
            Arc::new(chrom_arr),
            Arc::new(start_arr),
            Arc::new(strand_arr),
            Arc::new(best_cfd_arr),
            Arc::new(best_match_seq_arr),
        ],
    )?;

    let file = File::create(out_path)
        .with_context(|| format!("Failed to create output Parquet at {:?}", out_path))?;
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();

    let mut writer = ArrowWriter::try_new(file, schema, Some(props))?;
    writer.write(&batch)?;
    writer.close()?;

    Ok(())
}

/// Batch report (unchanged)
fn print_batch_report(
    batch_idx: usize,
    batch_rows: u64,
    read_ms: f64,
    stats: &BatchStats,
    running: &BatchStats,
    total_rows_so_far: u64,
    elapsed_search_s: f64,
) {
    let to_ms = |ns: u128| (ns as f64) / 1_000_000.0;

    let avg_bucket = stats.sum_bucket_depths.map(|s| (s as f64) / (batch_rows as f64));
    let avg_inter2 = (stats.sum_intersection2_size as f64) / (batch_rows as f64);
    let avg_inter3 = (stats.sum_intersection3_size as f64) / (batch_rows as f64);
    let avg_fallback_scan = (stats.sum_fallback_scanned as f64) / (batch_rows as f64);
    let avg_unique = (stats.sum_unique_scanned as f64) / (batch_rows as f64);
    let avg_hpass = (stats.sum_hamming_pass as f64) / (batch_rows as f64);
    let avg_cfd = (stats.sum_cfd_evals as f64) / (batch_rows as f64);
    let collect_ms = to_ms(stats.sum_collect_ns);
    let hamming_ms = to_ms(stats.sum_hamming_ns);

    let run_avg_bucket = running
        .sum_bucket_depths
        .map(|s| (s as f64) / (total_rows_so_far as f64));
    let run_avg_inter2 = (running.sum_intersection2_size as f64) / (total_rows_so_far as f64);
    let run_avg_inter3 = (running.sum_intersection3_size as f64) / (total_rows_so_far as f64);
    let run_avg_fallback = (running.sum_fallback_scanned as f64) / (total_rows_so_far as f64);
    let run_avg_unique = (running.sum_unique_scanned as f64) / (total_rows_so_far as f64);
    let run_avg_hpass = (running.sum_hamming_pass as f64) / (total_rows_so_far as f64);
    let run_avg_cfd = (running.sum_cfd_evals as f64) / (total_rows_so_far as f64);

    println!(
        "\n[Batch {batch_idx}] Rows: {batch_rows} | read_ms: {read_ms:.2} | collect_ms: {collect_ms:.2} | hamming_ms: {hamming_ms:.2}"
    );
    println!(
        "  Per-batch avg buckets (5): [{:.0}, {:.0}, {:.0}, {:.0}, {:.0}]",
        avg_bucket[0], avg_bucket[1], avg_bucket[2], avg_bucket[3], avg_bucket[4]
    );
    println!(
        "  Per-batch avg inter2: {:.1} | avg inter3: {:.1} | avg fallback_scanned: {:.1} | avg unique_scanned: {:.1}",
        avg_inter2, avg_inter3, avg_fallback_scan, avg_unique
    );
    println!(
        "  Per-batch avg hamming_pass<=r: {:.2} | avg CFD evals: {:.2}",
        avg_hpass, avg_cfd
    );
    println!(
        "  Running avg buckets:          [{:.0}, {:.0}, {:.0}, {:.0}, {:.0}]",
        run_avg_bucket[0], run_avg_bucket[1], run_avg_bucket[2], run_avg_bucket[3], run_avg_bucket[4]
    );
    println!(
        "  Running avg inter2: {:.1} | avg inter3: {:.1} | avg fallback_scanned: {:.1} | avg unique_scanned: {:.1} | elapsed: {:.1}s",
        run_avg_inter2, run_avg_inter3, run_avg_fallback, run_avg_unique, elapsed_search_s
    );
}

fn main() -> Result<()> {
    let args = Step3Args::parse();

    // Load CFD matrix from YAML (shared across workers)
    let cfd_matrix = load_cfd_matrix(&args.cfd_yaml, &args.cfd_name)?;
    let cfd_matrix = Arc::new(cfd_matrix);
    println!("CFD model loaded: '{}' (lsr_len={}) from {:?}",
        cfd_matrix.name, cfd_matrix.lsr_len, args.cfd_yaml
    );

    // Thread pool
    let available_cpus = rayon::current_num_threads();
    let active_threads = match args.threads {
        Some(t) if t > 0 => {
            let _ = rayon::ThreadPoolBuilder::new().num_threads(t).build_global();
            t
        }
        _ => available_cpus,
    };

    // Hardware info
    let mut sys = System::new_all();
    sys.refresh_all();
    let total_ram_bytes = sys.total_memory();

    println!("============================================================");
    println!("   GuideMaker Step-3 (Arrow, 3′) with YAML-loaded CFD + Parallel Indexing + 3-slice bitset   ");
    println!("============================================================");
    println!("Detected Hardware Configuration:");
    println!(" -> Active CPU Cores Available : {} Threads", available_cpus);
    println!(" -> System RAM Capacity        : {:.2} GB", total_ram_bytes as f64 / 1_073_741_824.0);
    println!("\nCalculated Optimization Plan:");
    println!(" -> Assigned Processing Threads: {}", active_threads);
    println!(" -> MIH Slice Params           : m = {} | L = {} nt | Offsets = {:?}", MIH_PARTITIONS, MIH_SLICE_LEN_NT, MIH_SLICE_OFFSETS);
    println!(" -> Prefilter Radius           : <= {} mismatches", MIH_MISMATCH_RADIUS);
    println!(" -> Arrow Batch Size           : {} rows/batch", args.batch_size);
    println!("============================================================");

    let start_pipeline = Instant::now();

    // Build candidate MIH index via Arrow streaming (3′ case)
    println!("\nBuilding candidate MIH index (candidate == true) via Arrow streaming...");
    let start_idx = Instant::now();
    let mut index = build_mih_index_arrow(&args.input, args.batch_size, args.max_bucket_depth)?;
    let index_build_time_sec = start_idx.elapsed().as_secs_f64();
    let total_candidates = index.candidates.len();

    // Bucket stats
    let total_entries: usize = index
        .buckets
        .iter()
        .map(|m| m.values().map(|b| b.ids.len()).sum::<usize>())
        .sum();
    let n_keys: usize = index.buckets.iter().map(|m| m.len()).sum();
    let avg_bucket_len = if n_keys > 0 {
        total_entries as f64 / n_keys as f64
    } else {
        0.0
    };
    let max_bucket_len = index
        .buckets
        .iter()
        .flat_map(|m| m.values().map(|b| b.ids.len()))
        .max()
        .unwrap_or(0);

    println!(
        "Candidate MIH Index Built: {} candidates indexed across {} bucket-keys in {:.2}s.",
        total_candidates, n_keys, index_build_time_sec
    );
    println!(
        "Bucket Statistics         : Avg Depth = {:.2} candidates/key | Max Depth = {} candidates",
        avg_bucket_len, max_bucket_len
    );

    // Stream targets via Arrow
    println!("\nStarting parallel streaming target search engine (Arrow, 3′)...");

    let file = File::open(&args.input)
        .with_context(|| format!("Failed to open input Parquet file at {:?}", args.input))?;
    let mut reader = ParquetRecordBatchReaderBuilder::try_new(file)?
        .with_batch_size(args.batch_size)
        .build()?; // Arrow RecordBatch reader

    let start_search_all = Instant::now();
    let mut batch_idx: usize = 0;
    let mut accum_search_time = 0.0;

    let mut running_stats = BatchStats::default();
    let mut total_rows_processed: u64 = 0;

    loop {
        let t_read_start = Instant::now();
        let maybe_batch = reader.next();
        let read_ms = t_read_start.elapsed().as_secs_f64() * 1000.0;

        let batch = match maybe_batch {
            Some(Ok(b)) => b,
            Some(Err(e)) => return Err(e.into()),
            None => break,
        };

        let start_chunk = Instant::now();
        let stats = process_targets_batch_against_index(
            &batch,
            &mut index,
            LSR_LEN,
            MIH_MISMATCH_RADIUS,
            &cfd_matrix,
        )?;
        let chunk_time = start_chunk.elapsed().as_secs_f64();
        accum_search_time += chunk_time;

        batch_idx += 1;
        total_rows_processed += stats.rows;

        running_stats.rows += stats.rows;
        for bi in 0..MIH_PARTITIONS {
            running_stats.sum_bucket_depths[bi] += stats.sum_bucket_depths[bi];
        }
        running_stats.sum_intersection2_size += stats.sum_intersection2_size;
        running_stats.sum_intersection3_size += stats.sum_intersection3_size;
        running_stats.sum_fallback_scanned += stats.sum_fallback_scanned;
        running_stats.sum_unique_scanned += stats.sum_unique_scanned;
        running_stats.sum_hamming_pass += stats.sum_hamming_pass;
        running_stats.sum_cfd_evals += stats.sum_cfd_evals;
        running_stats.sum_collect_ns += stats.sum_collect_ns;
        running_stats.sum_hamming_ns += stats.sum_hamming_ns;

        let elapsed_search = start_search_all.elapsed().as_secs_f64();

        // Per-batch report
        print_batch_report(
            batch_idx,
            stats.rows,
            read_ms,
            &stats,
            &running_stats,
            total_rows_processed,
            elapsed_search,
        );

        let global_rate = total_rows_processed as f64 / elapsed_search.max(1e-6);
        print!(
            "\r[Arrow Streaming] Batch {:>5} | Rows: {:>12} | Chunk Rate: {:>8.0} rows/s | Avg Rate: {:>8.0} rows/s | Elapsed: {:>6.1}s",
            batch_idx,
            total_rows_processed,
            (stats.rows as f64) / chunk_time.max(1e-6),
            global_rate,
            elapsed_search
        );
        let _ = std::io::stdout().flush();
    }

    println!(
        "\n\nProcessed {} Arrow batches ({} total targets) in {:.2}s.",
        batch_idx, total_rows_processed, accum_search_time
    );

    let probes_per_tgt = if total_rows_processed > 0 {
        index.total_bucket_probes as f64 / (total_rows_processed as f64)
    } else {
        0.0
    };
    let dedup_per_tgt = if total_rows_processed > 0 {
        index.total_candidates_deduped as f64 / (total_rows_processed as f64)
    } else {
        0.0
    };
    let cfd_per_cand = if total_candidates > 0 {
        index.total_cfd_evaluations as f64 / total_candidates as f64
    } else {
        0.0
    };

    println!("\n=== Search Engine Performance & Diagnostics Breakdown ===");
    println!(
        "Total Bucket Probes        : {} ({:.2} candidate probes/target)",
        index.total_bucket_probes, probes_per_tgt
    );
    println!(
        "Unique Candidates Checked  : {} ({:.2} deduped candidates/target)",
        index.total_candidates_deduped, dedup_per_tgt
    );
    println!(
        "CFD Matrix Evaluations     : {} ({:.2} CFD evaluations/candidate)",
        index.total_cfd_evaluations, cfd_per_cand
    );

    println!("\nWriting results to Parquet (candidates + best_* metrics)...");
    write_results_parquet(&args.out, &index)?;

    let global_total_time = start_pipeline.elapsed().as_secs_f64();
    println!("\n=== Step-3 (Arrow, 3′) Final Execution Summary ===");
    println!("Total Candidates Processed: {}", total_candidates);
    println!("Output saved successfully to: {:?}", args.out);
    println!(
        "Timings Breakdown         : Index Build: {:.2}s | Target Streaming: {:.2}s",
        index_build_time_sec, accum_search_time
    );
    println!("Total Walltime            : {:.2}s", global_total_time);
    println!("============================================================");

    Ok(())
}