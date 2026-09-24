use anyhow::{Context, Result};
use clap::Parser;
use guidemaker_scan::*;
use polars::prelude::*;
use std::fs::File;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "Step-3 Candidate Reduction via Off-Target Base Hamming Distance")]
pub struct Step3Args {
    /// Path to Step-2 guides Parquet file
    #[arg(long, alias = "input", visible_alias = "guides")]
    pub guides: PathBuf,

    /// Base Hamming distance threshold d (off-target mismatch limit). Candidates with a neighbor having distance < d are marked distpass=false.
    #[arg(long, short = 'd', default_value_t = 3)]
    pub d: u32,

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

    if let Some(num_threads) = args.threads {
        if num_threads > 0 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .ok();
        }
    }

    let guides_file = File::open(&args.guides)
        .with_context(|| format!("Failed to open guides Parquet file at {:?}", args.guides))?;
    let guides_df = ParquetReader::new(guides_file).finish()?;

    let (mut output_df, stats) = execute_step3(&guides_df, args.d, args.target_len)?;

    println!("=== Step-3 Hamming Distance Off-Target Summary ===");
    println!("Threshold d: {} base mismatches | Target Len: {} nt", args.d, args.target_len);
    println!("Total Candidate Inputs: {}", stats.total_candidates);
    println!("Passed Candidates (distpass=true): {}", stats.passed_candidates);
    println!("Failed Candidates (distpass=false): {}", stats.failed_candidates);
    println!("Pass Rate: {:.2}%", stats.pass_rate);
    println!(
        "Timings: Index Build: {:.4}s | MIH Distance Search: {:.4}s | Total Execution: {:.4}s",
        stats.index_build_time_sec, stats.search_time_sec, stats.total_time_sec
    );

    write_parquet(&mut output_df, &args.out)?;
    println!("Output written to {:?}", args.out);

    Ok(())
}
