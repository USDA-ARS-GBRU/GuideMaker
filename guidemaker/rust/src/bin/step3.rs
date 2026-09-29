use anyhow::{Context, Result};
use clap::Parser;
use guidemaker_scan::*;
use polars::prelude::*;
use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about = "Step-3 Candidate Reduction via Off-Target Base Hamming Distance & CFD Scoring")]
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

    /// Output Parquet file path
    #[arg(long)]
    pub out: PathBuf,
}

fn main() -> Result<()> {
    polars::enable_string_cache();
    let args = Step3Args::parse();

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

    println!("============================================================");
    println!("   GuideMaker Step-3 Full-Genome Off-Target Execution Engine ");
    println!("============================================================");
    println!("Detected Hardware Configuration:");
    println!(" -> Active CPU Cores Available : {} Threads", available_cpus);
    println!(" -> System RAM Capacity        : {:.2} GB", total_ram_bytes as f64 / 1_073_741_824.0);
    println!(" -> Assigned Processing Threads: {}", active_threads);
    println!("============================================================");

    let start_pipeline = Instant::now();

    let file = File::open(&args.guides)
        .with_context(|| format!("Failed to open input Parquet file at {:?}", args.guides))?;
    let guides_df = ParquetReader::new(file).finish()?;

    println!("Total Reference Targets Loaded : {}", guides_df.height());
    println!("Threshold Mismatch Limit     : d = {}", args.d);
    println!("Top-N Nearest Neighbors        : n = {}", args.top_n);

    let (mut final_df, stats) = execute_step3(
        &guides_df,
        args.d,
        args.top_n,
        &args.method,
        args.target_len,
    )?;

    let out_file = File::create(&args.out)
        .with_context(|| format!("Failed to create final output file at {:?}", args.out))?;
    ParquetWriter::new(out_file).finish(&mut final_df)?;

    let global_total_time = start_pipeline.elapsed().as_secs_f64();

    println!("\n=== Step-3 Final Summary ===");
    println!("Total Candidate Inputs : {}", stats.total_candidates);
    println!("Passed Candidates      : {}", stats.passed_candidates);
    println!("Failed Candidates      : {}", stats.failed_candidates);
    println!("Global Pass Rate       : {:.2}%", stats.pass_rate);
    println!(
        "Timings                : Index Build: {:.2}s | Nearest Neighbor Search: {:.2}s | Total: {:.2}s",
        stats.index_build_time_sec, stats.search_time_sec, stats.total_time_sec
    );
    println!("Total Pipeline Walltime: {:.2} Minutes", global_total_time / 60.0);
    println!("Pipeline Run Complete. Output saved to: {:?}", args.out);
    println!("============================================================");

    Ok(())
}
