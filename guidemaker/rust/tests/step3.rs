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
fn test_step3_slice_index_and_max_cfd() {
    let seq_cand1 = encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap();
    let seq_target_near = encode_2bit_u64(b"CCGTACGTACGTACGTACGT").unwrap(); // 1 mismatch (Pos 0: A->C)
    let seq_cand_nohit = encode_2bit_u64(b"TGCATGCATGCATGCATGCA").unwrap(); // 20 mismatches

    let hits = vec![
        TargetHit { candidate: true, seq: seq_cand1, chrom_idx: 0, start: 0, stop: 20, strand: true },
        TargetHit { candidate: false, seq: seq_target_near, chrom_idx: 0, start: 100, stop: 120, strand: true },
        TargetHit { candidate: true, seq: seq_cand_nohit, chrom_idx: 0, start: 200, stop: 220, strand: true },
    ];
    let chrom_names = vec!["chr1".to_string()];
    let df = build_dataframe(&hits, &chrom_names).unwrap();

    let (df_out, stats) = execute_step3(&df, 5, 3, "hnsw", 20).unwrap();
    assert_eq!(stats.total_candidates, 2);
    assert_eq!(df_out.height(), 2);

    assert!(df_out.get_column_names().iter().any(|&n| n == "best_target_seq_u64"));
    assert!(df_out.get_column_names().iter().any(|&n| n == "best_hamming"));
    assert!(df_out.get_column_names().iter().any(|&n| n == "best_cfd"));
    assert!(df_out.get_column_names().iter().any(|&n| n == "hits_scanned"));

    let seq_ca = df_out.column("seq").unwrap().u64().unwrap();
    let best_cfd_ca = df_out.column("best_cfd").unwrap().f32().unwrap();
    let best_ham_ca = df_out.column("best_hamming").unwrap().u32().unwrap();
    let best_tgt_ca = df_out.column("best_target_seq_u64").unwrap().u64().unwrap();
    let hits_ca = df_out.column("hits_scanned").unwrap().u32().unwrap();

    // Row 0: seq_cand1 -> should match seq_target_near
    assert_eq!(seq_ca.get(0), Some(seq_cand1));
    assert_eq!(best_tgt_ca.get(0), Some(seq_target_near));
    assert_eq!(best_ham_ca.get(0), Some(1));
    assert!(best_cfd_ca.get(0).unwrap() > 0.0);
    assert!(hits_ca.get(0).unwrap() >= 1);

    // Row 1: seq_cand_nohit -> no hit under prefilter threshold, should have sentinel defaults
    assert_eq!(seq_ca.get(1), Some(seq_cand_nohit));
    assert_eq!(best_tgt_ca.get(1), Some(0));
    assert_eq!(best_ham_ca.get(1), Some(255));
    assert_eq!(best_cfd_ca.get(1), Some(0.0));
    assert_eq!(hits_ca.get(1), Some(0));
}

#[test]
fn test_step3_integration_cli() {
    polars::enable_string_cache();
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("input_guides.parquet");
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
            candidate: false,
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
    write_parquet(&mut guides_df, &input_path).unwrap();

    let bin = env!("CARGO_BIN_EXE_guidemaker-step3");
    let status = Command::new(bin)
        .arg("--input")
        .arg(&input_path)
        .arg("--lsr-len")
        .arg("20")
        .arg("--slice-len")
        .arg("5")
        .arg("--slice-offsets")
        .arg("2,7,12")
        .arg("--prefilter-mismatch")
        .arg("5")
        .arg("--out")
        .arg(&out_step3)
        .status()
        .unwrap();

    assert!(status.success());
    assert!(out_step3.exists());

    let df_out = ParquetReader::new(File::open(&out_step3).unwrap())
        .finish()
        .unwrap();

    assert_eq!(df_out.height(), 2);
    assert!(df_out.get_column_names().iter().any(|name| name.as_str() == "best_target_seq_u64"));
    assert!(df_out.get_column_names().iter().any(|name| name.as_str() == "best_hamming"));
    assert!(df_out.get_column_names().iter().any(|name| name.as_str() == "best_cfd"));
    assert!(df_out.get_column_names().iter().any(|name| name.as_str() == "hits_scanned"));

    let cand_col = df_out.column("candidate").unwrap().bool().unwrap();
    assert_eq!(cand_col.get(0), Some(true));
    assert_eq!(cand_col.get(1), Some(true));
}
