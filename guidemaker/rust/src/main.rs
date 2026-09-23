use anyhow::{anyhow, Result};
use clap::Parser;
use guidemaker_scan::*;
use rayon::prelude::*;
use std::path::PathBuf;
use std::time::Instant;

/// CLI arguments for guidemaker-scan
#[derive(Parser, Debug)]
#[command(author, version, about = "CRISPR PAM & Target Scanner in Rust")]
pub struct Args {
    /// Path to FASTA or GenBank sequence file (plain, .gz, or .zst)
    #[arg(long, alias = "seq-file", visible_alias = "seq_file")]
    pub fasta: PathBuf,

    /// Optional path to GFF or GTF annotation file (plain, .gz, or .zst)
    #[arg(long, alias = "gtf", visible_alias = "annotation")]
    pub gff: Option<PathBuf>,

    /// PAM sequence (IUPAC ambiguous string, e.g. NGG)
    #[arg(long)]
    pub pam: String,

    /// Orientation: 5prime or 3prime
    #[arg(long)]
    pub orientation: String,

    /// Target length (1..=26, default 20)
    #[arg(long, default_value_t = 20)]
    pub target_len: usize,

    /// Number of worker threads (default: all available CPUs)
    #[arg(long)]
    pub threads: Option<usize>,

    /// Optional CSV output path for target hits
    #[arg(long)]
    pub out_csv: Option<PathBuf>,

    /// Optional Parquet output path for target hits
    #[arg(long)]
    pub out_parquet: Option<PathBuf>,

    /// Optional CSV output path for genomic features
    #[arg(long)]
    pub out_features_csv: Option<PathBuf>,

    /// Optional Parquet output path for genomic features
    #[arg(long)]
    pub out_features_parquet: Option<PathBuf>,
}

fn get_resource_usage() -> (f64, f64) {
    unsafe {
        let mut usage = std::mem::zeroed();
        if libc::getrusage(libc::RUSAGE_SELF, &mut usage) == 0 {
            let user_sec = usage.ru_utime.tv_sec as f64 + (usage.ru_utime.tv_usec as f64 / 1_000_000.0);
            let sys_sec = usage.ru_stime.tv_sec as f64 + (usage.ru_stime.tv_usec as f64 / 1_000_000.0);
            let cpu_sec = user_sec + sys_sec;

            #[cfg(target_os = "macos")]
            let rss_mb = (usage.ru_maxrss as f64) / (1024.0 * 1024.0);

            #[cfg(not(target_os = "macos"))]
            let rss_mb = (usage.ru_maxrss as f64) / 1024.0;

            (cpu_sec, rss_mb)
        } else {
            (0.0, 0.0)
        }
    }
}

fn main() -> Result<()> {
    let start_wall = Instant::now();

    let args = Args::parse();

    if let Some(num_threads) = args.threads {
        if num_threads > 0 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .ok();
        }
    }

    if !(1..=26).contains(&args.target_len) {
        return Err(anyhow!(
            "target-len must be between 1 and 26 inclusive, got {}",
            args.target_len
        ));
    }

    let orientation_is_5prime = match args.orientation.to_lowercase().as_str() {
        "5prime" => true,
        "3prime" => false,
        other => return Err(anyhow!("Invalid orientation: '{}'. Must be '5prime' or '3prime'", other)),
    };

    let pam_masks = parse_pam_masks(&args.pam)?;

    let records = read_sequence_records(&args.fasta)?;
    if records.len() > u16::MAX as usize {
        return Err(anyhow!("Sequence record count exceeds u16::MAX (65535)"));
    }

    let chrom_names: Vec<String> = records.iter().map(|r| r.id.clone()).collect();

    let all_hits: Vec<TargetHit> = records
        .par_iter()
        .enumerate()
        .flat_map(|(idx, record)| {
            let chrom_idx = idx as u32;
            let mut hits = Vec::new();

            if orientation_is_5prime {
                hits.extend(search_5prime_forward(chrom_idx, &record.seq, &pam_masks, args.target_len));
                hits.extend(search_5prime_reverse(chrom_idx, &record.seq, &pam_masks, args.target_len));
            } else {
                hits.extend(search_3prime_forward(chrom_idx, &record.seq, &pam_masks, args.target_len));
                hits.extend(search_3prime_reverse(chrom_idx, &record.seq, &pam_masks, args.target_len));
            }

            hits
        })
        .collect();

    let mut df = build_dataframe(&all_hits, &chrom_names)?;

    println!("Total target rows: {}", df.height());
    println!("{}", df.head(Some(5)));

    if let Some(csv_path) = &args.out_csv {
        write_csv(&mut df, csv_path)?;
    }

    if let Some(parquet_path) = &args.out_parquet {
        write_parquet(&mut df, parquet_path)?;
    }

    // Process Genomic Features if requested or available
    let feature_records = if let Some(gff_path) = &args.gff {
        read_feature_records(gff_path)?
    } else {
        // Try reading features directly from the sequence file (e.g. GenBank)
        read_feature_records(&args.fasta).unwrap_or_default()
    };

    if !feature_records.is_empty() {
        let mut feat_df = build_features_dataframe(&feature_records)?;
        println!("Total feature rows: {}", feat_df.height());
        println!("{}", feat_df.head(Some(5)));

        if let Some(feat_csv) = &args.out_features_csv {
            write_csv(&mut feat_df, feat_csv)?;
        }
        if let Some(feat_parquet) = &args.out_features_parquet {
            write_parquet(&mut feat_df, feat_parquet)?;
        }
    }

    let wall_sec = start_wall.elapsed().as_secs_f64();
    let (cpu_sec, peak_rss_mb) = get_resource_usage();

    if peak_rss_mb >= 1024.0 {
        println!(
            "Resource Usage: Wall time: {:.2}s | CPU time: {:.2}s | Peak RSS: {:.2} GB ({:.1} MB)",
            wall_sec,
            cpu_sec,
            peak_rss_mb / 1024.0,
            peak_rss_mb
        );
    } else {
        println!(
            "Resource Usage: Wall time: {:.2}s | CPU time: {:.2}s | Peak RSS: {:.2} MB",
            wall_sec,
            cpu_sec,
            peak_rss_mb
        );
    }

    Ok(())
}
