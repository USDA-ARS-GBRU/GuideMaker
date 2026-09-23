use anyhow::{Context, Result};
use clap::Parser;
use guidemaker_scan::*;
use polars::prelude::*;
use std::fs::File;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "Step-2 Candidate Reduction via Feature Windows + LSR")]
pub struct Step2Args {
    /// Path to Step-1 guides Parquet file
    #[arg(long)]
    pub guides: PathBuf,

    /// Path to Step-1 features Parquet file
    #[arg(long)]
    pub features: PathBuf,

    /// Upstream window in bp relative to feature TSS
    #[arg(long, default_value_t = 2000)]
    pub before: u32,

    /// Downstream window in bp into feature relative to TSS
    #[arg(long, default_value_t = 500)]
    pub into: u32,

    /// Seed length in nt used for LSR uniqueness test
    #[arg(long, default_value_t = 8)]
    pub lsr_len: usize,

    /// Output Parquet file path
    #[arg(long)]
    pub out: PathBuf,

    /// Optional comma-separated list of feature types to restrict spatial filtering (e.g. gene,transcript,exon)
    #[arg(long, value_delimiter = ',')]
    pub feature_types: Option<Vec<String>>,

    /// Try spatial proximity filter first (then LSR)
    #[arg(long, default_value_t = false)]
    pub fast_filter_first: bool,
}

fn main() -> Result<()> {
    let args = Step2Args::parse();

    let guides_file = File::open(&args.guides)
        .with_context(|| format!("Failed to open guides Parquet file at {:?}", args.guides))?;
    let guides_df = ParquetReader::new(guides_file).finish()?;

    let features_file = File::open(&args.features)
        .with_context(|| format!("Failed to open features Parquet file at {:?}", args.features))?;
    let features_df = ParquetReader::new(features_file).finish()?;

    let ftypes_ref = args.feature_types.as_deref();

    let (mut filtered_df, stats) = execute_step2(
        &guides_df,
        &features_df,
        args.before,
        args.into,
        args.lsr_len,
        ftypes_ref,
        args.fast_filter_first,
    )?;

    println!("=== Step-2 Filtering & Benchmark Summary ===");
    println!("Strategy: {}", if args.fast_filter_first { "Spatial First -> LSR" } else { "LSR First -> Spatial" });
    println!("Total Input Rows: {}", stats.total_input_rows);
    println!("Rows Passing Spatial Filter: {}", stats.rows_passing_spatial);
    println!("Rows Passing LSR Uniqueness: {}", stats.rows_passing_lsr);
    println!("Final Candidate Rows (Passing Both): {}", stats.final_candidate_rows);
    println!(
        "Timings: Spatial Filter: {:.4}s | LSR Uniqueness: {:.4}s | Total Execution: {:.4}s",
        stats.spatial_time_sec, stats.lsr_time_sec, stats.total_time_sec
    );

    write_parquet(&mut filtered_df, &args.out)?;
    println!("Output written to {:?}", args.out);

    Ok(())
}
