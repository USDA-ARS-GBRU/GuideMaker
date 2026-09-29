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

    // 1. Initialize modern sysinfo hardware tracking (SystemExt is no longer needed)
    let mut sys = sysinfo::System::new_all();
    sys.refresh_all();
    let total_ram_bytes = sys.total_memory(); 

    // 2. Query Rayon or the host hardware layout to determine available CPU cores natively
    let available_cpus = rayon::current_num_threads();
    let active_threads = match args.threads {
        Some(t) => {
            if t > 0 {
                // Reinitialize the global thread pool to match user constraints if explicitly passed
                let _ = rayon::ThreadPoolBuilder::new().num_threads(t).build_global();
            }
            t
        }
        None => available_cpus,
    };

    // 3. AUTOMATED CHUNK TUNING:
    // Check if the user left chunk_size at its default value (10_000_000). 
    // If they didn't pass a custom override, dynamically calculate an optimal bound.
    if args.chunk_size == 10_000_000 {
        // Approximate bytes consumed per 1 million active rows given top_n matrix layers
        let bytes_per_million_rows = 50_000_000 * (args.top_n as u64);
        
        // Target an isolated safe ceiling per chunk pass (5% of total host RAM)
        let target_chunk_memory_allocation = total_ram_bytes / 20;

        // Calculate rows, bounding it tightly between 2M and 25M elements
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
    
    // ... The rest of your streaming file processing loop continues down below unchanged



    println!("============================================================");
    println!("   GuideMaker Step-3 Managed Memory Streaming Pipeline       ");
    println!("============================================================");

    let start_pipeline = Instant::now();

    // 1. Inspect Parquet metadata to calculate streaming boundaries without reading data
    let file = File::open(&args.guides)
        .with_context(|| format!("Failed to open input Parquet file at {:?}", args.guides))?;
    let mut reader = ParquetReader::new(file);
    let total_rows = reader.num_rows().context("Failed to read Parquet row metadata")?;
    
    println!("Total Rows Detected on Disk : {}", total_rows);
    println!("Streaming Chunk Size        : {}", args.chunk_size);
    println!("Threshold Target Mismatches : d = {}", args.d);

    let mut current_row = 0;
    let mut chunk_idx = 0;
    let mut temporary_file_paths = Vec::new();
    
    // Track global counters across individual chunk passes
    let mut accum_total_candidates = 0;
    let mut accum_passed_candidates = 0;
    let mut accum_failed_candidates = 0;
    let mut accum_idx_build_time = 0.0;
    let mut accum_search_time = 0.0;

    // 2. Core Slicing Loop: Isolate chunks sequentially out of scope to clear RAM
    while current_row < total_rows {
        let n_rows_to_read = std::cmp::min(args.chunk_size, total_rows - current_row);
        println!(
            "-> Processing Chunk {} | Rows {}-{}...",
            chunk_idx, current_row, current_row + n_rows_to_read
        );

        // Leverage the Lazy engine to scan and slice the file chunk cleanly from disk without memory spikes
        let chunk_df = LazyFrame::scan_parquet(&args.guides, ScanArgsParquet::default())?
            .slice(current_row as i64, n_rows_to_read as u32)
            .collect()?;

        // Execute our highly optimized zero-allocation Pigeonhole + Nearest Neighbor extraction
        let (mut processed_chunk_df, stats) = execute_step3(
            &chunk_df,
            args.d,
            args.top_n,
            &args.method,
            args.target_len,
        )?;

        // Update aggregation telemetry counters
        accum_total_candidates += stats.total_candidates;
        accum_passed_candidates += stats.passed_candidates;
        accum_failed_candidates += stats.failed_candidates;
        accum_idx_build_time += stats.index_build_time_sec;
        accum_search_time += stats.search_time_sec;

        // Write intermediate slice frame to disk to clear memory registers
        let tmp_path = format!("final_chunk_{}.parquet", chunk_idx);
        let tmp_file = File::create(&tmp_path)?;
        ParquetWriter::new(tmp_file).finish(&mut processed_chunk_df)?;

        temporary_file_paths.push(tmp_path);
        current_row += n_rows_to_read;
        chunk_idx += 1;
    }

    // 3. Assemble Master Dataset via Lazy Construction without loading data into memory
    println!("\nAll chunks processed. Assembling final unified output dataset via Lazy Schema...");
    
    let mut chunk_lazy_frames = Vec::new();
    for path in &temporary_file_paths {
        let lf = LazyFrame::scan_parquet(path, ScanArgsParquet::default())?;
        chunk_lazy_frames.push(lf);
    }

    // Vertically concatenate the lazy frame streams
    let joined_lazy = concat(chunk_lazy_frames, UnionArgs::default())?;
    
    let mut final_output_df = joined_lazy.collect().context("Failed to collect final concatenated dataframes")?;
    
    // Write out the pristine master output file
    let out_file = File::create(&args.out)
        .with_context(|| format!("Failed to create final output file at {:?}", args.out))?;
    ParquetWriter::new(out_file).finish(&mut final_output_df)?;

    // 4. Garbage collection: Erase intermediate scratch chunk files from the storage drive
    for path in temporary_file_paths {
        let _ = std::fs::remove_file(path);
    }

    let global_total_time = start_pipeline.elapsed().as_secs_f64();
    let global_pass_rate = if accum_total_candidates > 0 {
        (accum_passed_candidates as f64 / accum_total_candidates as f64) * 100.0
    } else {
        0.0
    };

    println!("\n=== Step-3 Final Consolidated Streaming Summary ===");
    println!("Total Candidate Inputs : {}", accum_total_candidates);
    println!("Passed Candidates      : {}", accum_passed_candidates);
    println!("Failed Candidates      : {}", accum_failed_candidates);
    println!("Global Pass Rate       : {:.2}%", global_pass_rate);
    println!(
        "Accumulated Timings    : Index Build: {:.2}s | Nearest Neighbor Search: {:.2}s",
        accum_idx_build_time, accum_search_time
    );
    println!("Total Pipeline Walltime: {:.2} Minutes", global_total_time / 60.0);
    println!("Pipeline Run Complete. Output saved to: {:?}", args.out);
    println!("============================================================");

    Ok(())
}
