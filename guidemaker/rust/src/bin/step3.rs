use anyhow::{Context, Result};
use clap::Parser;
use guidemaker_scan::*;
use polars::prelude::*;
use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about = "Step-3 Candidate Reduction via Off-Target Base Hamming Distance (Memory-Managed Streaming)")]
pub struct Step3Args {
    /// Path to Step-2 guides Parquet file
    #[arg(long, alias = "input", visible_alias = "guides")]
    pub guides: PathBuf,

    /// Base Hamming distance threshold d (off-target mismatch limit).
    #[arg(long, short = 'd', default_value_t = 3)]
    pub d: u32,

    /// Number of top nearest neighbors to compute for passing candidates
    #[arg(long, short = 'n', default_value_t = 3)]
    pub top_n: usize,

    /// Neighbor search method ("hnsw" or "exact")
    #[arg(long, default_value = "hnsw")]
    pub method: String,

    /// Target length in nt (1..=26, default 20)
    #[arg(long, default_value_t = 20)]
    pub target_len: usize,

    /// Number of worker threads (default: all available CPUs)
    #[arg(long)]
    pub threads: Option<usize>,

    /// Streaming batch chunk size (number of rows to process at once)
    #[arg(long, default_value_t = 10_000_000)]
    pub chunk_size: usize,

    /// Output Parquet file path
    #[arg(long)]
    pub out: PathBuf,
}

fn main() -> Result<()> {
    polars::enable_string_cache();
    let mut args = Step3Args::parse();

    let mut sys = sysinfo::System::new_all();
    sys.refresh_all();
    let total_ram_bytes = sys.total_memory();

    let available_cpus = rayon::current_num_threads();
    let active_threads = match args.threads {
        Some(t) => {
            if t > 0 {
                let _ = rayon::ThreadPoolBuilder::new().num_threads(t).build_global();
            }
            t
        }
        None => available_cpus,
    };

    if args.chunk_size == 10_000_000 {
        let bytes_per_million_rows = 50_000_000 * (args.top_n as u64);
        let target_chunk_memory_allocation = total_ram_bytes / 20;
        let calculated_rows = (target_chunk_memory_allocation / bytes_per_million_rows) * 1_000_000;
        args.chunk_size = calculated_rows.clamp(2_000_000, 25_000_000) as usize;
    }

    println!("============================================================");
    println!("   GuideMaker Step-3 Auto-Tuning Execution Engine           ");
    println!("============================================================");
    println!("Detected Hardware Configuration:");
    println!(" -> Active CPU Cores Available : {} Threads", available_cpus);
    println!(" -> System RAM Capacity        : {:.2} GB", total_ram_bytes as f64 / 1_073_741_824.0);
    println!("\nCalculated Optimization Plan:");
    println!(" -> Assigned Processing Threads: {}", active_threads);
    println!(" -> Automated Chunk Boundary   : {} rows per batch", args.chunk_size);
    println!("============================================================");

    println!("============================================================");
    println!("   GuideMaker Step-3 Managed Memory Streaming Pipeline       ");
    println!("============================================================");

    let start_pipeline = Instant::now();

    let file = File::open(&args.guides)
        .with_context(|| format!("Failed to open input Parquet file at {:?}", args.guides))?;
    let mut reader = ParquetReader::new(file);
    let total_rows = reader.num_rows().context("Failed to read Parquet row metadata")?;

    println!("Total Reference Targets Loaded : {}", total_rows);
    println!("Streaming Chunk Size           : {}", args.chunk_size);
    println!("Threshold Mismatch Limit       : d = {}", args.d);
    println!("Top-N Nearest Neighbors        : n = {}", args.top_n);

    // 1. Scan and load ONLY candidate guide rows (candidate == true) into memory (~100MB RAM)
    println!("\nLoading Candidate Guides (candidate == true)...");
    let cand_df = LazyFrame::scan_parquet(&args.guides, ScanArgsParquet::default())?
        .filter(col("candidate").eq(lit(true)))
        .collect()?;

    println!("Building Candidate Guide Index for {} candidate guides...", cand_df.height());
    let start_idx = Instant::now();
    let candidate_index = build_candidate_index(&cand_df, args.top_n)?;
    let index_build_time_sec = start_idx.elapsed().as_secs_f64();
    let total_candidates = candidate_index.candidates.len();
    println!(
        "Candidate Index Built: {} candidate guides indexed in {:.2}s.",
        total_candidates, index_build_time_sec
    );

    // 2. Stream all genomic targets in auto-tuned chunks from disk against candidate index
    let mut current_row = 0;
    let mut chunk_idx = 0;
    let mut accum_search_time = 0.0;

    while current_row < total_rows {
        let n_rows_to_read = std::cmp::min(args.chunk_size, total_rows - current_row);
        println!(
            "-> Scanning Genomic Targets Chunk {} | Rows {}-{}...",
            chunk_idx, current_row, current_row + n_rows_to_read
        );

        let targets_chunk_df = LazyFrame::scan_parquet(&args.guides, ScanArgsParquet::default())?
            .slice(current_row as i64, n_rows_to_read as u32)
            .collect()?;

        let start_chunk = Instant::now();
        scan_targets_against_candidates(&targets_chunk_df, &candidate_index, args.d, args.target_len)?;
        accum_search_time += start_chunk.elapsed().as_secs_f64();

        current_row += n_rows_to_read;
        chunk_idx += 1;
    }

    // 3. Compile final results DataFrame for passing candidate guides
    println!("\nCompiling final output dataset for passing candidate guides...");
    let mut final_df = build_candidate_results_dataframe(&cand_df, candidate_index, args.top_n, args.target_len)?;

    let passed_candidates = final_df.height();
    let failed_candidates = total_candidates.saturating_sub(passed_candidates);
    let pass_rate = if total_candidates > 0 {
        (passed_candidates as f64 / total_candidates as f64) * 100.0
    } else {
        0.0
    };

    let out_file = File::create(&args.out)
        .with_context(|| format!("Failed to create final output file at {:?}", args.out))?;
    ParquetWriter::new(out_file).finish(&mut final_df)?;

    let global_total_time = start_pipeline.elapsed().as_secs_f64();

    println!("\n=== Step-3 Final Consolidated Streaming Summary ===");
    println!("Total Candidate Inputs : {}", total_candidates);
    println!("Passed Candidates      : {}", passed_candidates);
    println!("Failed Candidates      : {}", failed_candidates);
    println!("Global Pass Rate       : {:.2}%", pass_rate);
    println!(
        "Accumulated Timings    : Index Build: {:.2}s | Nearest Neighbor Search: {:.2}s",
        index_build_time_sec, accum_search_time
    );
    println!("Total Pipeline Walltime: {:.2} Minutes", global_total_time / 60.0);
    println!("Pipeline Run Complete. Output saved to: {:?}", args.out);
    println!("============================================================");

    Ok(())
}
