use anyhow::{Context, Result};
use std::io::Read;

/// Representation of a single parsed sequence record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqRecord {
    pub id: String,
    pub seq: Vec<u8>,
}

/// Parse FASTA records from reader
pub fn parse_fasta_records<R: Read>(reader: R) -> Result<Vec<SeqRecord>> {
    let fasta_reader = bio::io::fasta::Reader::new(reader);
    let mut records = Vec::new();
    for rec_res in fasta_reader.records() {
        let rec = rec_res.context("Failed to parse FASTA record")?;
        records.push(SeqRecord {
            id: rec.id().to_string(),
            seq: rec.seq().to_ascii_uppercase(),
        });
    }
    Ok(records)
}

/// Parse GenBank records from reader
pub fn parse_genbank_records<R: Read>(reader: R) -> Result<Vec<SeqRecord>> {
    let buf_reader = std::io::BufReader::new(reader);
    let mut records = Vec::new();

    let mut current_id = String::new();
    let mut current_seq = Vec::new();
    let mut in_origin = false;

    for line_res in buf_reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();

        if trimmed.starts_with("LOCUS") {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() > 1 {
                current_id = parts[1].to_string();
            } else {
                current_id = "unknown_contig".to_string();
            }
            current_seq.clear();
            in_origin = false;
        } else if trimmed.starts_with("ACCESSION") && current_id.is_empty() {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() > 1 {
                current_id = parts[1].to_string();
            }
        } else if trimmed.starts_with("VERSION") && current_id.is_empty() {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() > 1 {
                current_id = parts[1].to_string();
            }
        } else if trimmed.starts_with("ORIGIN") {
            in_origin = true;
        } else if trimmed == "//" {
            if !current_id.is_empty() && !current_seq.is_empty() {
                records.push(SeqRecord {
                    id: std::mem::take(&mut current_id),
                    seq: std::mem::take(&mut current_seq),
                });
            } else {
                current_id.clear();
                current_seq.clear();
            }
            in_origin = false;
        } else if in_origin {
            for b in line.bytes() {
                if b.is_ascii_alphabetic() {
                    current_seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }

    if !current_id.is_empty() && !current_seq.is_empty() {
        records.push(SeqRecord {
            id: current_id,
            seq: current_seq,
        });
    }

    Ok(records)
}

/// Read FASTA or GenBank sequence records transparently
pub fn read_sequence_records(path: &std::path::Path) -> Result<Vec<SeqRecord>> {
    let mut reader = super::open_compressed_reader(path)?;

    let mut header_buf = [0u8; 1024];
    let n = reader.read(&mut header_buf)?;

    let peek_str = String::from_utf8_lossy(&header_buf[..n]);

    let combined_reader = std::io::Cursor::new(header_buf[..n].to_vec()).chain(reader);

    if peek_str.trim_start().starts_with('>') {
        parse_fasta_records(combined_reader)
    } else if peek_str.contains("LOCUS") || peek_str.contains("ORIGIN") || peek_str.contains("FEATURES") {
        parse_genbank_records(combined_reader)
    } else {
        parse_fasta_records(combined_reader)
    }
}
