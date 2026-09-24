use anyhow::Result;
use clap::Parser;
use guidemaker_scan::{execute_step3, build_dataframe, TargetHit, encode_2bit_u64};
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(name = "benchmark-human-scale")]
#[command(about = "Generate synthetic human-scale guide dataset and profile Step-3 MIH performance")]
struct Args {
    /// Number of synthetic candidate target guides to generate
    #[arg(short = 'n', long, default_value_t = 1_000_000)]
    num_guides: usize,

    /// Off-target Hamming distance threshold d
    #[arg(short = 'd', default_value_t = 2)]
    d: u32,

    /// Target guide length in nt
    #[arg(long = "target-len", default_value_t = 20)]
    target_len: usize,

    /// Number of Rayon threads
    #[arg(short = 't', long, default_value_t = 8)]
    threads: usize,
}

fn generate_synthetic_dna(idx: usize, target_len: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seq = Vec::with_capacity(target_len);
    let mut state = idx as u64;
    for _ in 0..target_len {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let base = bases[(state >> 62) as usize];
        seq.push(base);
    }
    seq
}

fn main() -> Result<()> {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    println!("============================================================");
    println!("   GuideMaker Step-3 Benchmark Utility (Human-Scale Synthetic)");
    println!("============================================================");
    println!("Generating {} synthetic target candidates (target_len={} nt)...", args.num_guides, args.target_len);

    let start_gen = Instant::now();
    let hits: Vec<TargetHit> = (0..args.num_guides)
        .map(|i| {
            let dna = generate_synthetic_dna(i, args.target_len);
            let encoded = encode_2bit_u64(&dna).unwrap();
            TargetHit {
                candidate: true,
                seq: encoded,
                chrom_idx: (i % 24) as u32,
                start: (i * 50) as u32,
                stop: (i * 50 + args.target_len) as u32,
                strand: i % 2 == 0,
            }
        })
        .collect();

    let chrom_names: Vec<String> = (0..24).map(|i| format!("chr{}", i + 1)).collect();
    let df = build_dataframe(&hits, &chrom_names)?;
    println!("Generated DataFrame in {:.3}s", start_gen.elapsed().as_secs_f64());

    println!("\nExecuting Step-3 MIH Off-Target Hamming Filter (d={}, threads={})...", args.d, args.threads);
    let (_out_df, stats) = execute_step3(&df, args.d, args.target_len)?;

    println!("\n--- Step-3 Benchmark Results ---");
    println!("Total Candidate Guides : {}", stats.total_candidates);
    println!("Passed Candidates      : {}", stats.passed_candidates);
    println!("Failed Candidates      : {}", stats.failed_candidates);
    println!("Pass Rate              : {:.2}%", stats.pass_rate);
    println!("MIH Index Build Time   : {:.4} s", stats.index_build_time_sec);
    println!("MIH Search Engine Time : {:.4} s", stats.search_time_sec);
    println!("Total Step-3 Time      : {:.4} s", stats.total_time_sec);
    println!("Throughput             : {:.0} guides/sec", stats.total_candidates as f64 / stats.total_time_sec);
    println!("============================================================");

    Ok(())
}
