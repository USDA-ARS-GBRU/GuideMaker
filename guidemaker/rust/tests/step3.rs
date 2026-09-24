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
fn test_8_2_segment_partition_bounds() {
    // For L=20, d=2 -> m = d + 1 = 3 segments
    // L=20 / 3 = 6, rem = 2 -> segments bounds:
    // s=0: (0, 7)
    // s=1: (7, 7)
    // s=2: (14, 6)
    let m3 = 3;
    let (s0_start, s0_len) = get_segment_bounds(0, m3, 20);
    let (s1_start, s1_len) = get_segment_bounds(1, m3, 20);
    let (s2_start, s2_len) = get_segment_bounds(2, m3, 20);

    assert_eq!((s0_start, s0_len), (0, 7));
    assert_eq!((s1_start, s1_len), (7, 7));
    assert_eq!((s2_start, s2_len), (14, 6));
    assert_eq!(s0_len + s1_len + s2_len, 20);

    // For L=25, d=3 -> m = d + 1 = 4 segments
    // L=25 / 4 = 6, rem = 1 -> segments bounds:
    // s=0: (0, 7)
    // s=1: (7, 6)
    // s=2: (13, 6)
    // s=3: (19, 6)
    let m4 = 4;
    let (_s0_25_start, s0_25_len) = get_segment_bounds(0, m4, 25);
    let (_s1_25_start, s1_25_len) = get_segment_bounds(1, m4, 25);
    let (_s2_25_start, s2_25_len) = get_segment_bounds(2, m4, 25);
    let (_s3_25_start, s3_25_len) = get_segment_bounds(3, m4, 25);

    assert_eq!(s0_25_len + s1_25_len + s2_25_len + s3_25_len, 25);
    assert_eq!([s0_25_len, s1_25_len, s2_25_len, s3_25_len], [7, 6, 6, 6]);
}

#[test]
fn test_8_3_mih_candidate_evaluation() {
    let mask_20 = compute_target_mask(20);

    // Synthetic dataset of sequences
    let mut seq_vec = Vec::new();
    let seq_base = encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap();
    seq_vec.push(seq_base); // ID 0

    // Sequence with 1 mismatch to seq_base (pos 0: A vs C)
    let seq_near1 = encode_2bit_u64(b"CCGTACGTACGTACGTACGT").unwrap();
    seq_vec.push(seq_near1); // ID 1

    // Sequence with 5 mismatches to seq_base
    let seq_far = encode_2bit_u64(b"TGCATGCATGCATGCATGCA").unwrap();
    seq_vec.push(seq_far); // ID 2

    // Build MIH Index for d=2 (m=3)
    let index = MihIndex::build(&seq_vec, 20, 2);

    // Candidate 0 vs seq_near1 (distance 1 < d=2) -> should fail
    let pass_0 = index.evaluate_candidate(0, seq_base, &seq_vec, 2, mask_20);
    assert_eq!(pass_0, false);

    // Candidate 2 vs seq_base and seq_near1 (distances >= 5 > d=2) -> should pass
    let pass_2 = index.evaluate_candidate(2, seq_far, &seq_vec, 2, mask_20);
    assert_eq!(pass_2, true);
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

    assert_eq!(df_out.height(), 3);
    assert!(df_out.get_column_names().iter().any(|name| name.as_str() == "distpass"));

    let distpass_ca = df_out.column("distpass").unwrap().bool().unwrap();
    assert_eq!(distpass_ca.get(0), Some(false));
    assert_eq!(distpass_ca.get(1), Some(false));
    assert_eq!(distpass_ca.get(2), Some(true));
}
