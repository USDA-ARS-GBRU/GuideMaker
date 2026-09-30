# Specification: GuideMaker Step-3 High-Throughput Search Pipeline

## 1. Objective & Core Task
Rebuild the Step-3 candidate reduction and Top-N nearest-neighbor search pipeline. The engine must 
evaluate a small pool of high-value candidate guides against a massive streaming set 
of genomic target sites to remove candidates celow a thhreshold hamming distance and for the remaining canidates, identify the N nearst neighbors and their hammin distances

## 2. Key Data Dimensions & Constraints
* **Query Set (Candidates):** ~1,645,409 rows (loaded  as a polars dataframe with 2-bit encoded uint64 ).
* **Reference Database (Targets):** ~317,350,744 rows (streamed sequentially from disk in batches of ~8,000,000 rows dunamically set based on resources).
* **DNA sequence** the originating data is a DNA string {G,A,T,C} with a target length of 20-27 characters, plus ad pam sit of 3 plus two additional charaters, total lenth does not exceed 32. this is 2 bit encoded. The additioan data is passed for doench efficiency estimation but is masked for  hamming distance opperations and fitlering
* **Hardware Profile:** 15 active CPU threads, 24.00 GB System RAM ceiling.
* **Target Metric:** The final pipeline must complete execution over all 317M targets in minutes (e.g., 1 to 60 min), faster is better.

## 3. Mandatory Architectural Criteria

### A. Inverted Index Mapping (Candidates Indexed, Targets Streamed)
* **Rule:** Do NOT index the 317M target database. It breaks the RAM ceiling and causes massive performance degradation.
* **Implementation:** Build the Pigeonhole Index exclusively on the 1.6M candidate guides extracted from teh polars dataframe bsased on the candidate flag. This keeps the index size under ~50MB, pinning it tightly within the CPU's L3/L2 cache.
* **Streaming Loop:** Process the 317M targets by passing sequential chunks (e.g., 8M rows at a time) into the scanning engine.

### B. Mathematical Parity in Bit Extraction (Pigeonhole Constraint)
* **Rule:** The block extraction logic must map precisely to the target window length without discarding bits or causing bucket collisions.
* **20nt Spacer Logic (Bits 24..63):** Use 4 independent 10-bit blocks (5 bases each) ( modify number of bins and block size to optimize search efficenty for other lengths. you may need to determinethis experimentally):
  ```rust
  pub fn extract_pigeonhole_blocks(seq: u64) -> [usize; 4] {
      [
          ((seq >> 54) & 0x03FF) as usize, // Block 0
          ((seq >> 44) & 0x03FF) as usize, // Block 1
          ((seq >> 34) & 0x03FF) as usize, // Block 2
          ((seq >> 24) & 0x03FF) as usize, // Block 3
      ]
  }
  ```
* **Distance Boundary:** A 4-block setup mathematically guarantees finding targets where Hamming distance \(d \le 3\). If your application configuration sets \(d > 3\), the number of blocks must be scaled dynamically to \(d + 1\) chunks to prevent false negatives.

### C. Zero-Allocation Parallel Execution Loop
* **Rule:** The parallel iterator (`into_par_iter()`) running across targets must be completely free of heap allocations (`Vec::new()`, `.collect()`, cloning strings) to eliminate global allocator lock contention.
* **Grain Size:** Configure Rayon chunks using `.with_min_len(4096)` to minimize thread context-switching overhead.

### D. Contention-Free Thread Synchronization
* **Fast-Failing:** Candidate records must use an lock-free atomic gate (`AtomicBool`) to track viability:
  ```rust
  if !cand.passed.load(Ordering::Relaxed) { return; }
  ```
* **Lock Elimination:** Do NOT call a heavy OS mutex lock for standard target iterations. Use `AtomicBool` for fast termination when an off-target with \(dist < d\) is encountered. For valid nearest-neighbor collection (\(dist \ge d\)), protect the collection array with a non-blocking spinlock or a trylock routine (`try_lock()`).

### E. Positional Row & Vector Alignment
* **Rule:** When compiling the final Polars DataFrame, results vectors (`nn_dist`, `nn_seq`, `nn_cfd`) must strictly correspond to their original row alignment inside the input polars dataframe.
* **Implementation:** Map results directly back via an immutable tracker index (`guide_df_idx`) to avoid vector shifting or row duplication bugs.

## 4. Evaluation & Performance Metrics
Before merging, verify the codebase against the following performance checkpoints:
1. **CPU Scaling:** Core utilization across your 15 threads should remain near 100% during the target scanning phase (no threads hanging or sleeping due to lock contention).
2. **Memory Footprint:** Resident memory usage must remain flat and predictable beneath the 24.00 GB ceiling throughout all 40 streaming chunks.
3. **Deterministic Alignment:** Confirm that running the pipeline yields identical rows and correctly structured list series outputs where metrics perfectly match their target candidate rows.


### 5. Output 

* The script should retun a polars dataframe in parquet format  with the candidates after candidates with an hamming distance below threshold have been removed. the data frame should have new columns: 
1. neighbors: a list with the top N targets in 2bit encoded uint64 format.
2. nn_dist: a list with the hamming distance for the targets in the same order as the neighbor list.
3. cfd: A list with the CDF scores of the nearist neighbors.