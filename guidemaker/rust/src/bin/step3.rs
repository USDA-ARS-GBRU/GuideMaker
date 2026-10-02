use anyhow::{Context, Result};
use clap::Parser;
use guidemaker_scan::*;
use polars::prelude::*;
use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about = "Step-3 Candidate Reduction via Slice-Based Inverted Index & Max-CFD Off-Target Search")]
pub struct Step3Args {
    /// Path to input Parquet file (unified target/guide dataset)
    #[arg(long, visible_alias = "guides", visible_alias = "input")]
    pub input: PathBuf,

    /// Deprecated threshold / prefilter mismatch limit (alias for --prefilter-mismatch)
    #[arg(long, short = 'd')]
    pub d: Option<u32>,

    /// Deprecated parameter (ignored)
    #[arg(long, short = 'n')]
    pub top_n: Option<usize>,

    /// Deprecated parameter (ignored)
    #[arg(long)]
    pub method: Option<String>,

    /// Deprecated parameter (alias for --lsr-len)
    #[arg(long)]
    pub target_len: Option<usize>,

    /// LSR length for Hamming filter & CFD evaluation (default 20 nt)
    #[arg(long, default_value_t = 20)]
    pub lsr_len: usize,

    /// Length L of inverted index key slices in nt (default 5)
    #[arg(long, default_value_t = 5)]
    pub slice_len: usize,

    /// Comma-separated 0-based start offset positions for slices within LSR (default "2,7,12")
    #[arg(long, default_value = "2,7,12")]
    pub slice_offsets: String,

    /// Maximum LSR mismatches allowed for prefilter evaluation (default 5)
    #[arg(long, default_value_t = 5)]
    pub prefilter_mismatch: u32,

    /// Number of worker threads (default: all available CPUs)
    #[arg(long)]
    pub threads: Option<usize>,

    /// Streaming batch chunk size (number of target rows per chunk)
    #[arg(long, default_value_t = 100_000)]
    pub batch_size: usize,

    /// Whether to keep non-candidate rows in output DataFrame (default false)
    #[arg(long, default_value_t = false)]
    pub keep_noncandidates: bool,

    /// Output Parquet file path
    #[arg(long)]
    pub out: PathBuf,
}

fn main() -> Result<()> {
    polars::enable_string_cache();
    let args = Step3Args::parse();

    let slice_offsets_vec: Vec<usize> = args
        .slice_offsets
        .split(',')
        .map(|s| s.trim().parse::<usize>().context("Failed to parse slice_offsets"))
        .collect::<Result<Vec<usize>>>()?;

    let effective_lsr_len = args.target_len.unwrap_or(args.lsr_len);
    let effective_prefilter_mismatch = args.d.unwrap_or(args.prefilter_mismatch);

    let config = Step3Config {
        lsr_len: effective_lsr_len,
        slice_len: args.slice_len,
        slice_offsets: slice_offsets_vec,
        prefilter_mismatch: effective_prefilter_mismatch,
        keep_noncandidates: args.keep_noncandidates,
    };
    config.validate()?;

    let mut sys = sysinfo::System::new_all();
    sys.refresh_all();
    let total_ram_bytes = sys.total_memory();

    let available_cpus = rayon::current_num_threads();
    let active_threads = match args.threads {
        Some(t) if t > 0 => {
            let _ = rayon::ThreadPoolBuilder::new().num_threads(t).build_global();
            t
        }
        _ => available_cpus,
    };

    println!("============================================================");
    println!("   GuideMaker Step-3 Slice Index & Max-CFD Execution Engine ");
    println!("============================================================");
    println!("Detected Hardware Configuration:");
    println!(" -> Active CPU Cores Available : {} Threads", available_cpus);
    println!(" -> System RAM Capacity        : {:.2} GB", total_ram_bytes as f64 / 1_073_741_824.0);
    println!("\nCalculated Optimization Plan:");
    println!(" -> Assigned Processing Threads: {}", active_threads);
    println!(" -> Inverted Index Slice Params : L = {} nt | Offsets = {:?}", config.slice_len, config.slice_offsets);
    println!(" -> Prefilter Mismatch Threshold: < {} mismatches", config.prefilter_mismatch);
    println!(" -> Target Streaming Batch Size: {} rows per batch", args.batch_size);
    println!("============================================================");

    let start_pipeline = Instant::now();

    let file = File::open(&args.input)
        .with_context(|| format!("Failed to open input Parquet file at {:?}", args.input))?;
    let mut reader = ParquetReader::new(file);
    let total_rows = reader.num_rows().context("Failed to read Parquet row metadata")?;

    println!("Total Reference Targets Loaded : {}", total_rows);

    // 1. Scan input table and extract candidate guides (candidate == true)
    println!("\nScanning input dataset for candidate guides (candidate == true)...");
    let full_df = LazyFrame::scan_parquet(&args.input, ScanArgsParquet::default())?.collect()?;

    let start_idx = Instant::now();
    let index = build_step3_candidate_index(&full_df, &config)?;
    let index_build_time_sec = start_idx.elapsed().as_secs_f64();
    let total_candidates = index.candidates.len();
    println!(
        "Candidate Slice Index Built: {} candidate guides indexed in {:.2}s.",
        total_candidates, index_build_time_sec
    );

    // 2. Stream all target rows in batch chunks against candidate slice index
    let mut current_row = 0;
    let mut chunk_idx = 0;
    let mut accum_search_time = 0.0;

    while current_row < total_rows {
        let n_rows_to_read = std::cmp::min(args.batch_size, total_rows - current_row);
        let targets_chunk_df = LazyFrame::scan_parquet(&args.input, ScanArgsParquet::default())?
            .slice(current_row as i64, n_rows_to_read as u32)
            .collect()?;

        let start_chunk = Instant::now();
        scan_targets_batch_against_index(&targets_chunk_df, &index, &config)?;
        accum_search_time += start_chunk.elapsed().as_secs_f64();

        current_row += n_rows_to_read;
        chunk_idx += 1;
    }

    println!("\nProcessed {} streaming target batches in {:.2}s.", chunk_idx, accum_search_time);

    // 3. Materialize final DataFrame with appended best_* columns
    println!("Compiling final output dataset with max-CFD off-target metrics...");
    let mut final_df = build_step3_results_dataframe(&full_df, index, &config)?;

    let out_file = File::create(&args.out)
        .with_context(|| format!("Failed to create final output file at {:?}", args.out))?;
    ParquetWriter::new(out_file).finish(&mut final_df)?;

    let global_total_time = start_pipeline.elapsed().as_secs_f64();

    println!("\n=== Step-3 Final Execution Summary ===");
    println!("Total Candidates Processed: {}", total_candidates);
    println!("Output Rows Written       : {}", final_df.height());
    println!(
        "Timings Breakdown         : Index Build: {:.2}s | Target Streaming: {:.2}s",
        index_build_time_sec, accum_search_time
    );
    println!("Total Walltime            : {:.2}s", global_total_time);
    println!("Output saved successfully to: {:?}", args.out);
    println!("============================================================");

    Ok(())
}
