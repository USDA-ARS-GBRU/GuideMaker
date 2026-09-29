use guidemaker_scan::*;
use polars::prelude::*;
use std::fs::File;
use std::process::Command;
use tempfile::tempdir;

#[test]
fn test_8_1_base_hamming_correctness() {
    let mask_20 = compute_target_mask(20);

    // Identical sequences
    let seq_a = encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap();
    let seq_b = encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap();
    assert_eq!(base_hamming_distance_masked(seq_a, seq_b, mask_20), 0);

    // Mismatches at specific positions
    // Pos 0: A vs C
    let seq_pos0 = encode_2bit_u64(b"CCGTACGTACGTACGTACGT").unwrap();
    assert_eq!(base_hamming_distance_masked(seq_a, seq_pos0, mask_20), 1);

    // Pos 19 (last base): T vs G
    let seq_pos19 = encode_2bit_u64(b"ACGTACGTACGTACGTACGG").unwrap();
    assert_eq!(base_hamming_distance_masked(seq_a, seq_pos19, mask_20), 1);

    // 3 bases changed (pos 0, 1, 2)
    let seq_3diff = encode_2bit_u64(b"TTTTACGTACGTACGTACGT").unwrap();
    assert_eq!(base_hamming_distance_masked(seq_a, seq_3diff, mask_20), 3);

    // Fully different sequence (all 20 bases different)
    let seq_all_diff = encode_2bit_u64(b"TGCATGCATGCATGCATGCA").unwrap();
    assert_eq!(base_hamming_distance_masked(seq_a, seq_all_diff, mask_20), 20);
}

#[test]
fn test_8_3_step3_hnsw_and_exact_methods() {
    let seq_base = encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap();
    let seq_near = encode_2bit_u64(b"CCGTACGTACGTACGTACGT").unwrap(); // 1 mismatch
    let seq_far = encode_2bit_u64(b"TGCATGCATGCATGCATGCA").unwrap(); // 20 mismatches

    let hits = vec![
        TargetHit { candidate: true, seq: seq_base, chrom_idx: 0, start: 0, stop: 20, strand: true },
        TargetHit { candidate: true, seq: seq_near, chrom_idx: 0, start: 100, stop: 120, strand: true },
        TargetHit { candidate: true, seq: seq_far, chrom_idx: 0, start: 200, stop: 220, strand: true },
    ];
    let chrom_names = vec!["chr1".to_string()];
    let df = build_dataframe(&hits, &chrom_names).unwrap();

    // Test Exact method
    let (df_exact, stats_exact) = execute_step3(&df, 2, 2, "exact", 20).unwrap();
    assert_eq!(stats_exact.passed_candidates, 1);
    assert_eq!(df_exact.height(), 1);
    assert_eq!(df_exact.column("seq").unwrap().u64().unwrap().get(0), Some(seq_far));
    assert!(df_exact.get_column_names().iter().any(|&n| n == "nn_dist"));
    assert!(df_exact.get_column_names().iter().any(|&n| n == "nn_seq"));

    // Test HNSW method
    let (df_hnsw, stats_hnsw) = execute_step3(&df, 2, 2, "hnsw", 20).unwrap();
    assert_eq!(stats_hnsw.passed_candidates, 1);
    assert_eq!(df_hnsw.height(), 1);
    assert_eq!(df_hnsw.column("seq").unwrap().u64().unwrap().get(0), Some(seq_far));
    assert!(df_hnsw.get_column_names().iter().any(|&n| n == "nn_dist"));
    assert!(df_hnsw.get_column_names().iter().any(|&n| n == "nn_seq"));
}

#[test]
fn test_8_4_step3_integration_cli() {
    polars::enable_string_cache();
    let dir = tempdir().unwrap();
    let guides_path = dir.path().join("guides_input.parquet");
    let out_step3 = dir.path().join("out_step3.parquet");

    let hits = vec![
        TargetHit {
            candidate: true,
            seq: encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap(),
            chrom_idx: 0,
            start: 1000,
            stop: 1020,
            strand: true,
        },
        TargetHit {
            candidate: true,
            seq: encode_2bit_u64(b"ACGTACGTACGTACGTACGA").unwrap(), // 1 mismatch to hit 0
            chrom_idx: 0,
            start: 2000,
            stop: 2020,
            strand: true,
        },
        TargetHit {
            candidate: true,
            seq: encode_2bit_u64(b"TGCATGCATGCATGCATGCA").unwrap(), // 10 mismatches to hits 0 and 1
            chrom_idx: 0,
            start: 3000,
            stop: 3020,
            strand: true,
        },
    ];
    let chrom_names = vec!["chr1".to_string()];
    let mut guides_df = build_dataframe(&hits, &chrom_names).unwrap();
    write_parquet(&mut guides_df, &guides_path).unwrap();

    let bin = env!("CARGO_BIN_EXE_guidemaker-step3");
    let status = Command::new(bin)
        .arg("--guides")
        .arg(&guides_path)
        .arg("-d")
        .arg("2")
        .arg("-n")
        .arg("2")
        .arg("--method")
        .arg("hnsw")
        .arg("--target-len")
        .arg("20")
        .arg("--out")
        .arg(&out_step3)
        .status()
        .unwrap();

    assert!(status.success());
    assert!(out_step3.exists());

    let df_out = ParquetReader::new(File::open(&out_step3).unwrap())
        .finish()
        .unwrap();

    assert_eq!(df_out.height(), 1);
    assert!(df_out.get_column_names().iter().any(|name| name.as_str() == "nn_dist"));
    assert!(df_out.get_column_names().iter().any(|name| name.as_str() == "nn_seq"));
    assert_eq!(df_out.column("candidate").unwrap().bool().unwrap().get(0), Some(true));
}
