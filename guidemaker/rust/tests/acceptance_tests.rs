use guidemaker_scan::*;
use polars::prelude::*;
use std::fs::File;
use std::io::Write;
use std::process::Command;
use std::sync::Arc;
use tempfile::tempdir;

#[test]
fn acceptance_test_1_iupac_and_orientation() {
    let seq = b"ACGTACGTACGTACGTACGTCGG";
    let pam_masks = parse_pam_masks("NGG").unwrap();
    let chrom_names = vec!["chr1".to_string()];

    let fwd_hits = search_3prime_forward(0, seq, &pam_masks, 20);
    assert_eq!(fwd_hits.len(), 1);
    assert_eq!(fwd_hits[0].candidate, true);
    assert_eq!(fwd_hits[0].chrom_idx, 0);
    assert_eq!(fwd_hits[0].start, 0);
    assert_eq!(fwd_hits[0].stop, 20);
    assert_eq!(fwd_hits[0].strand, true);
    assert_eq!(fwd_hits[0].orientation, false);

    let df = build_dataframe(&fwd_hits, &chrom_names).unwrap();
    assert_eq!(df.height(), 1);

    let expected_code = encode_2bit_u64(b"ACGTACGTACGTACGTACGT").unwrap();
    assert_eq!(fwd_hits[0].seq, expected_code);
}

#[test]
fn acceptance_test_2_reverse_mapping_correctness() {
    let fwd_seq = b"CCTACGTACGTACGTACGTACGT";
    let pam_masks = parse_pam_masks("NGG").unwrap();

    let rev_hits = search_3prime_reverse(0, fwd_seq, &pam_masks, 20);
    assert_eq!(rev_hits.len(), 1);
    assert_eq!(rev_hits[0].candidate, true);
    assert_eq!(rev_hits[0].chrom_idx, 0);
    assert_eq!(rev_hits[0].start, 3);
    assert_eq!(rev_hits[0].stop, 23);
    assert_eq!(rev_hits[0].strand, false);
    assert_eq!(rev_hits[0].orientation, false);

    assert!(rev_hits[0].start < rev_hits[0].stop);
    assert!(rev_hits[0].stop <= fwd_seq.len() as u32);
}

#[test]
fn acceptance_test_3_encoding_check() {
    let seq = b"ACGTACGTACGTACGTACGT";
    let encoded = encode_2bit_u64(seq).unwrap();
    assert_eq!(encoded, 0x1B1B_1B1B_1B00_0000);
}

#[test]
fn acceptance_test_4_csv_parquet_parity_and_cli() {
    let dir = tempdir().unwrap();
    let fasta_path = dir.path().join("input.fasta");
    let csv_path = dir.path().join("out.csv");
    let parquet_path = dir.path().join("out.parquet");

    let mut f = File::create(&fasta_path).unwrap();
    writeln!(f, ">chr0").unwrap();
    writeln!(f, "ACGTACGTACGTACGTACGTCGGTTTTACGTACGTACGTACGTACGTCGG").unwrap();

    let bin = env!("CARGO_BIN_EXE_guidemaker-scan");
    let status = Command::new(bin)
        .arg("--fasta")
        .arg(&fasta_path)
        .arg("--pam")
        .arg("NGG")
        .arg("--orientation")
        .arg("3prime")
        .arg("--target-len")
        .arg("20")
        .arg("--out-csv")
        .arg(&csv_path)
        .arg("--out-parquet")
        .arg(&parquet_path)
        .status()
        .unwrap();

    assert!(status.success());
    assert!(csv_path.exists());
    assert!(parquet_path.exists());

    let parquet_df = ParquetReader::new(File::open(&parquet_path).unwrap())
        .finish()
        .unwrap();

    let mut schema = Schema::default();
    schema.insert_at_index(0, "candidate".into(), DataType::Boolean).unwrap();
    schema.insert_at_index(1, "seq".into(), DataType::UInt64).unwrap();
    schema.insert_at_index(2, "chrom".into(), DataType::Categorical(None, Default::default())).unwrap();
    schema.insert_at_index(3, "start".into(), DataType::UInt32).unwrap();
    schema.insert_at_index(4, "stop".into(), DataType::UInt32).unwrap();
    schema.insert_at_index(5, "orientation".into(), DataType::Boolean).unwrap();
    schema.insert_at_index(6, "strand".into(), DataType::Boolean).unwrap();

    let csv_df = CsvReadOptions::default()
        .with_schema(Some(Arc::new(schema)))
        .try_into_reader_with_file_path(Some(csv_path))
        .unwrap()
        .finish()
        .unwrap();

    assert_eq!(csv_df.height(), parquet_df.height());
    assert_eq!(csv_df.width(), parquet_df.width());
    assert_eq!(csv_df.get_column_names(), parquet_df.get_column_names());
    assert_eq!(csv_df.schema(), parquet_df.schema());
}

#[test]
fn acceptance_test_5_boundary_conditions_and_ambiguous_skipping() {
    let short_seq = b"CGG";
    let pam_masks = parse_pam_masks("NGG").unwrap();
    let hits = search_3prime_forward(0, short_seq, &pam_masks, 20);
    assert!(hits.is_empty());

    let seq_with_n = b"ACGTACGTACGTACGTACNTCGG";
    let hits_n = search_3prime_forward(0, seq_with_n, &pam_masks, 20);
    assert!(hits_n.is_empty(), "Target containing 'N' should be skipped");
}
