# GuideMaker Step‑3 — Jules Replacement Spec (Concise)

## 1) Objective

Write a standa alone rust script that evaluates each **candidate guide** (from a Polars DataFrame, 2‑bit encoded a subset of all rows) against a **streaming set of genomic targets** to:
1. filter the candidates-target pairs to those with a hamming distance of 5 or less and
2. for each sett passign hte filter score the pair  using CDF scoring and  retain the maximum score and the sequence of the max scoring  target for each cadndidate.

**Core rules**

* Use a **slice‑based inverted index** + **bitwise Hamming filter** to minimize CFD work.
* Hamming filter: **20‑nt LSR**, **mismatches < 5** pass; otherwise skip CFD.
* For passing pairs: compute **CFD (Doench)**; track **highest CFD** per candidate with the associated **2‑bit sequence** of the target.
* Designed to stream \~317M targets quickly with \~15 threads and ≤24 GB RAM.

***

## 2) Inputs & Data Model

** data frame example**

```
   candidate                   seq         chrom  start   stop  strand
0      False   1567635393809849024  NC_000001.11  10445  10470    True
1      False   6270548258482969168  NC_000001.11  10458  10483    True
2      False  12418475214316853936  NC_000001.11  10471  10496    True
3      False  12780412709848312528  NC_000001.11  10472  10497    True
4      False  11918043498572909088  NC_000001.11  10484  10509    True
```

* Seq is encoded as [target] [pam] [pad] of 3prime  orientation, [pad] [pam] [target] if 5prime oreentation. LSR is the area closest to the pam on the target. 
* the LSR is the first 20 nt nearest the pam site withing the target.
* Each candidate is unique in hte 20 nt lsr  but some targets may be exact duplicated
 

**Candidates (queries)**

* \~1.9M guides; where the candidates column is True in a Polars DataFrame.
* 2‑bit encoded sequences in `u64` (A=00, C=01, G=10, T=11).
* LSR length for Hamming filter: **20 nt** (default).

**Targets (reference)**

* \~317M target sequences, 2‑bit encoded `u64`.
* Streamed in batches (from Parquet/Arrow/flat binary).

**CFD model**

* Source data `guidemaker/data/cfd_data.json` → compact LUT `[f32; 320]` (positions × base pairs).  you must precompute this once, now into fast binary lookup format.
* CFD score = product of per‑mismatch penalties.

**Hardware profile**

* Threads: \~15; RAM budget: ≤24 GB.

***

## 3) Outputs

A Polars DataFrame (Parquet) with **one row per candidate** (where `candidate=true`) containing:

* `best_target_seq_u64` — the 2‑bit encoded **target** that produced the highest CFD.
* `best_hamming` — Hamming distance (0–4) on the 20‑nt LSR.
* `best_cfd` — highest CFD score observed.
* `hits_scanned` — number of candidate–target pairs that passed the prefilter.

*Note:* This replaces prior “top‑N neighbors”; here we retain **only the best (max CFD)** per candidate.

***

## 4) High‑Level Algorithm

1. **Build slice‑based inverted index (Humti‑style)**
   * Choose **slice length L = 5** and **K = 3 slices** at fixed offsets in the 20‑nt LSR (e.g., positions 2, 7, 12).
   * For each candidate:
     * Extract the K slices → direct‑address keys in `[0, 2^(2L))`.
     * Append candidate ID to each slice’s bucket.
   * Result: `bucket_index[NUM_BUCKETS]` + contiguous payload of candidate IDs (and optionally sequences).

2. **Stream targets and probe buckets**
   * For each target:
     * Compute the same K slice keys from its LSR.
     * Collect the **union** of candidate IDs from matching buckets (dedup).

3. **Bitwise prefilter (Hamming on 20 nt)**
   * For each candidate in the union:
     * `x = cand_u64 ^ tgt_u64`
     * `mismatches = popcount(x) / 2`
     * If `mismatches >= 5` → **skip**.

4. **CFD scoring (mismatches only)**
   * Walk mismatch lanes only (see §6) and multiply LUT penalties:
     * `CFD = ∏ weight(pos, q_base, t_base)`
   * If `CFD > current_best[cand]` → **update**:
     * `best_target_seq_u64`, `best_hamming`, `best_cfd`.

5. **Write output**
   * After full pass, materialize the best results per candidate to Parquet.

***

## 5) Index Details (Direct‑Address Buckets)

* **Parameters**: `L=5` → `NUM_BUCKETS = 1 << (2*L) = 1024`.
* **Two‑pass build**:
  1. Count entries per bucket for all candidate slices.
  2. Prefix‑sum to set bucket starts; fill a single contiguous payload (IDs, optionally sequences).
* **Query**:
  * For a target, union candidate IDs from `K` buckets; dedup (sort‑dedup or thread‑local set).
* **Why it works**:
  * If two sequences differ at <5 positions, they typically share ≥1 contiguous window (L≥4), so the bucket probe prunes most candidates.

***

## 6) Bitwise & CFD Details

**2‑bit Hamming**

```text
x = candidate ^ target           // XOR
mismatches = popcount(x) / 2     // each base uses 2 bits
if mismatches >= 5 -> skip CFD
```

**CFD lookup indexing**

* For mismatch at position `pos` (0-based within 20‑nt LSR):
  * `q_base = (candidate >> (2*pos)) & 0b11`
  * `t_base = (target   >> (2*pos)) & 0b11`
  * `idx = (pos << 4) + (q_base << 2) + t_base`
  * `weight = CFD_LUT[idx]` (zero/missing ⇒ treat as 1.0)
* **Score**: multiply weights over mismatch positions only.

**Mismatch‑only walker (branch‑light)**

* Use `tzcnt`/`trailing_zeros` on `x` to find the next set two‑bit lane, clear it, and continue until `x==0` or `pos>=LSR_LEN`.

***

## 7) Minimal API / CLI (pick one)


**Optional CLI**

```text
guidemaker-step3-jules \
  --guides <parquet> \
  --targets <parquet|bin> \
  --lsr-len 20 --slice-len 5 --slice-offsets 2,7,12 \
  --prefilter-mismatch 5 --threads 15 --batch-size 100000 \
  --out step3_best.parquet
```

***

## 8) Efficient SIMD‑Friendly / Multithreaded Implementation (Short Guide)

1. **Direct‑address buckets**
   * Keys: `0..(1<<(2*L))-1`.
   * Use `Vec<Bucket>` + contiguous payload vectors (IDs and optionally candidate sequences).

2. **Contiguous candidate loads**
   * Duplicate candidate **sequences** per bucket (not just IDs) to enable linear scans:
     * `x = cand ^ target` → `popcount` → mismatch count.
   * Hardware prefetchers love contiguous reads; this reduces random memory probes.

3. **Branch‑free mismatch walker**
   * Skip 0..LSR\_LEN loops; iterate only set lanes:
     * `tzcnt`/`trailing_zeros` → lane index
     * clear lane bits → next
   * Multiply LUT penalties using precomputed base indices.

4. **Precomputed tables**
   * `SHIFT[pos] = 2*pos`, `IDX_BASE[pos] = pos<<4`.
   * Tiny static arrays accelerate base extraction and LUT addressing.

5. **Thread‑local scratch**
   * Keep per‑thread union buffers and neighbor stores; reuse memory to avoid allocator churn.

6. **Work‑stealing parallelism**
   * Chunk targets into **50k–150k** sub‑batches and feed a work‑stealing pool (e.g., Rayon/Crossbeam in Rust).
   * Evens out load from variable bucket sizes; keeps cores saturated.

7. **u64 everywhere**
   * Stick to XOR, shifts, masks, POPCNT/TZCNT; these map cleanly to SIMD/CPU instructions.

**Tuning defaults**

* `L=5`, `K=3`, `prefilter_mismatch=5`, `batch_size=100k`, `threads=15`.

***

## 9) Correctness & Validation

* **Unit tests**
  * 2‑bit encoding/decoding (A,C,G,T mapping).
  * Hamming on the 20‑nt LSR (XOR+popcount).
  * CFD indexing and score multiplication (spot‑check known pairs).

* **Integration**
  * Small synthetic dataset: ensure best target per candidate matches ground truth (within tolerance).
  * Determinism: fixed seeds, stable slice offsets, same output across runs.

* **Sanity checks**
  * Mismatch distribution histograms post‑prefilter (expect vast majority skipped).
  * CFD ranges and outliers; log candidate IDs with unusually high counts of near neighbors.

***

## 10) Resource & Edge‑Case Notes

* **Memory**: duplicating sequences per bucket costs \~`C * K * 8 bytes` (C=candidates, K=slices).  
  For `C≈1.9M`, `K=3` → \~43 MB for sequences (+ID payload, index metadata) — acceptable within 24 GB.

* **Sequence length**: focus Hamming on **exactly the 20‑nt LSR**. If guides/targets are longer (PAM/padding), do not include those bases in the mismatch filter; CFD penalties should use the same 20‑nt window unless model specifies otherwise.

* **Missing LUT entries**: treat as `1.0` (neutral) to avoid penalizing unsupported pairs.

* **Deduping unions**: use sort‑dedup for small unions; fall back to a thread‑local set for larger unions.

***

## 11) Minimal Pseudocode (hot path)

```rust
for target_batch in stream_targets(batch_size) {
  parallel_for target in target_batch {
    cand_ids = union_buckets(slice_keys(target, L, OFFSETS));  // dedup
    for cand_id in cand_ids {
      cand_seq = cand_seq_u64[cand_id];
      let x = cand_seq ^ target;
      let mismatches = x.count_ones() >> 1;
      if mismatches >= 5 { continue; }

      // mismatch-only CFD
      let mut y = x;
      let mut cfd = 1.0f32;
      while y != 0 {
        let tz = y.trailing_zeros() as usize;
        let pos = tz >> 1;
        if pos >= LSR_LEN { break; }
        y &= !(0b11u64 << (2*pos));
        let q = ((cand_seq >> (2*pos)) & 0b11) as usize;
        let t = ((target   >> (2*pos)) & 0b11) as usize;
        let w = CFD_LUT[(pos<<4) + (q<<2) + t];
        cfd *= if w == 0.0 { 1.0 } else { w };
      }
      if cfd > best_cfd[cand_id] {
        best_cfd[cand_id] = cfd;
        best_target_seq[cand_id] = target;
        best_hamming[cand_id] = mismatches as u8;
      }
    }
  }
}
write_parquet(guides_df.with_columns(best_*), out_path);
```

***

## 12) Deliverables

1. Script or CLI implementing the above with:
   * Index build (two‑pass buckets).
   * Streaming target scanning.
   * XOR+popcount prefilter (`<5` mismatches).
   * CFD mismatch‑only walker.
   * Best‑per‑candidate tracking.
   * Polars output as Parquet.

2. README (1‑pager):
   * How to run; parameter defaults; performance tips.

*