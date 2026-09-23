use anyhow::{anyhow, Context, Result};
use bio::io::fasta;
use clap::Parser;
use guidemaker_scan::*;
use std::path::PathBuf;

/// CLI arguments for guidemaker-scan
#[derive(Parser, Debug)]
#[command(author, version, about = "CRISPR PAM & Target Scanner in Rust")]
pub struct Args {
    /// Path to FASTA file
    #[arg(long)]
    pub fasta: PathBuf,

    /// PAM sequence (IUPAC ambiguous string, e.g. NGG)
    #[arg(long)]
    pub pam: String,

    /// Orientation: 5prime or 3prime
    #[arg(long)]
    pub orientation: String,

    /// Target length (1..=26, default 20)
    #[arg(long, default_value_t = 20)]
    pub target_len: usize,

    /// Optional CSV output path
    #[arg(long)]
    pub out_csv: Option<PathBuf>,

    /// Optional Parquet output path
    #[arg(long)]
    pub out_parquet: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();

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

    let reader = fasta::Reader::from_file(&args.fasta)
        .with_context(|| format!("Failed to open FASTA file at {:?}", args.fasta))?;

    let mut all_hits = Vec::new();
    let mut chrom_count: u64 = 0;

    for record_res in reader.records() {
        let record = record_res.with_context(|| "Failed to read FASTA record")?;
        if chrom_count > u16::MAX as u64 {
            return Err(anyhow!("FASTA record count exceeds u16::MAX (65535)"));
        }
        let chrom_idx = chrom_count as u16;
        chrom_count += 1;

        let seq_uppercase = record.seq().to_ascii_uppercase();

        if orientation_is_5prime {
            all_hits.extend(search_5prime_forward(chrom_idx, &seq_uppercase, &pam_masks, args.target_len));
            all_hits.extend(search_5prime_reverse(chrom_idx, &seq_uppercase, &pam_masks, args.target_len));
        } else {
            all_hits.extend(search_3prime_forward(chrom_idx, &seq_uppercase, &pam_masks, args.target_len));
            all_hits.extend(search_3prime_reverse(chrom_idx, &seq_uppercase, &pam_masks, args.target_len));
        }
    }

    let mut df = build_dataframe(&all_hits)?;

    println!("Total rows: {}", df.height());
    println!("{}", df.head(Some(5)));

    if let Some(csv_path) = &args.out_csv {
        write_csv(&mut df, csv_path)?;
    }

    if let Some(parquet_path) = &args.out_parquet {
        write_parquet(&mut df, parquet_path)?;
    }

    Ok(())
}
