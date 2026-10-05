use anyhow::{anyhow, Context, Result};
use bio::alphabets::dna;
use polars::prelude::*;

/// Representation of a single candidate hit
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TargetHit {
    pub candidate: bool,
    pub seq: u64,
    pub chrom_idx: u32,
    pub start: u32,
    pub stop: u32,
    pub strand: bool,
}

/// Convert an IUPAC character to its bitmask.
/// A=1, C=2, G=4, T=8.
pub fn iupac_mask(c: char) -> Result<u8> {
    match c.to_ascii_uppercase() {
        'A' => Ok(1),
        'C' => Ok(2),
        'G' => Ok(4),
        'T' => Ok(8),
        'M' => Ok(1 | 2),
        'R' => Ok(1 | 4),
        'W' => Ok(1 | 8),
        'S' => Ok(2 | 4),
        'Y' => Ok(2 | 8),
        'K' => Ok(4 | 8),
        'V' => Ok(1 | 2 | 4),
        'H' => Ok(1 | 2 | 8),
        'D' => Ok(1 | 4 | 8),
        'B' => Ok(2 | 4 | 8),
        'X' | 'N' => Ok(1 | 2 | 4 | 8),
        other => Err(anyhow!("Invalid IUPAC character: '{}'", other)),
    }
}

/// Convert a PAM string into a vector of per-base IUPAC bitmasks.
pub fn parse_pam_masks(pam: &str) -> Result<Vec<u8>> {
    if pam.is_empty() {
        return Err(anyhow!("PAM string cannot be empty"));
    }
    pam.chars().map(iupac_mask).collect()
}

/// Check if PAM matches the sequence at position `pos`.
pub fn pam_matches_forward(pam_masks: &[u8], seq_bytes: &[u8], pos: usize) -> bool {
    if pos + pam_masks.len() > seq_bytes.len() {
        return false;
    }
    for (i, &pam_mask) in pam_masks.iter().enumerate() {
        let seq_char = seq_bytes[pos + i] as char;
        let seq_mask = match iupac_mask(seq_char) {
            Ok(m) => m,
            Err(_) => return false,
        };
        if (pam_mask & seq_mask) == 0 {
            return false;
        }
    }
    true
}

/// Encode a target sequence (up to 26 nt) into a left-aligned 2-bit u64.
/// A=00, C=01, G=10, T=11.
/// Returns error if non-ATCG base encountered or length outside 1..=26.
pub fn encode_2bit_u64(seq_bytes: &[u8]) -> Result<u64> {
    if seq_bytes.is_empty() || seq_bytes.len() > 26 {
        return Err(anyhow!(
            "Target length must be between 1 and 26, got {}",
            seq_bytes.len()
        ));
    }
    let mut val: u64 = 0;
    for (i, &b) in seq_bytes.iter().enumerate() {
        let code: u64 = match b.to_ascii_uppercase() {
            b'A' => 0b00,
            b'C' => 0b01,
            b'G' => 0b10,
            b'T' => 0b11,
            other => return Err(anyhow!("Non-ATCG base encountered in target: '{}'", other as char)),
        };
        let shift = 64 - 2 * (i + 1);
        val |= code << shift;
    }
    Ok(val)
}

/// Reverse complement wrapper using rust-bio
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    dna::revcomp(seq)
}

/// Search 5prime orientation on forward strand: 5'-[PAM][TARGET]-3'
pub fn search_5prime_forward(
    chrom_idx: u32,
    seq: &[u8],
    pam_masks: &[u8],
    target_len: usize,
) -> Vec<TargetHit> {
    let mut hits = Vec::new();
    let pam_len = pam_masks.len();
    let flank_len = 2;
    let seq_len = seq.len();
    let total_window = target_len + pam_len + flank_len;
    if seq_len < total_window {
        return hits;
    }
    let min_pos = flank_len;
    let max_pos = seq_len - pam_len - target_len;
    for pos in min_pos..=max_pos {
        if pam_matches_forward(pam_masks, seq, pos) {
            let window_start = pos - flank_len;
            let window_end = pos + pam_len + target_len;
            let window_bytes = &seq[window_start..window_end];
            if let Ok(encoded_seq) = encode_2bit_u64(window_bytes) {
                let target_start = pos + pam_len;
                let target_end = target_start + target_len;
                let start = target_start as u32;
                let stop = target_end as u32;
                if start < stop && (stop as usize) <= seq_len {
                    hits.push(TargetHit {
                        candidate: true,
                        seq: encoded_seq,
                        chrom_idx,
                        start,
                        stop,
                        strand: true,
                    });
                }
            }
        }
    }
    hits
}

/// Search 5prime orientation on reverse strand
pub fn search_5prime_reverse(
    chrom_idx: u32,
    seq: &[u8],
    pam_masks: &[u8],
    target_len: usize,
) -> Vec<TargetHit> {
    let mut hits = Vec::new();
    let pam_len = pam_masks.len();
    let flank_len = 2;
    let seq_len = seq.len();
    let total_window = target_len + pam_len + flank_len;
    if seq_len < total_window {
        return hits;
    }
    let rc_seq = dna::revcomp(seq);
    let min_rc_pos = flank_len;
    let max_rc_pos = seq_len - pam_len - target_len;
    for rc_pos in min_rc_pos..=max_rc_pos {
        if pam_matches_forward(pam_masks, &rc_seq, rc_pos) {
            let rc_window_start = rc_pos - flank_len;
            let rc_window_end = rc_pos + pam_len + target_len;
            let window_bytes = &rc_seq[rc_window_start..rc_window_end];
            if let Ok(encoded_seq) = encode_2bit_u64(window_bytes) {
                let start_fwd = seq_len - (rc_pos + pam_len + target_len);
                let stop_fwd = seq_len - (rc_pos + pam_len);
                if start_fwd < stop_fwd && stop_fwd <= seq_len {
                    hits.push(TargetHit {
                        candidate: true,
                        seq: encoded_seq,
                        chrom_idx,
                        start: start_fwd as u32,
                        stop: stop_fwd as u32,
                        strand: false,
                    });
                }
            }
        }
    }
    hits
}

/// Search 3prime orientation on forward strand: 5'-[TARGET][PAM]-3'
pub fn search_3prime_forward(
    chrom_idx: u32,
    seq: &[u8],
    pam_masks: &[u8],
    target_len: usize,
) -> Vec<TargetHit> {
    let mut hits = Vec::new();
    let pam_len = pam_masks.len();
    let flank_len = 2;
    let seq_len = seq.len();
    let total_window = target_len + pam_len + flank_len;
    if seq_len < total_window {
        return hits;
    }
    let min_pos = target_len;
    let max_pos = seq_len - pam_len - flank_len;
    for pos in min_pos..=max_pos {
        if pam_matches_forward(pam_masks, seq, pos) {
            let target_start = pos - target_len;
            let window_end = target_start + total_window;
            let window_bytes = &seq[target_start..window_end];
            if let Ok(encoded_seq) = encode_2bit_u64(window_bytes) {
                let start = target_start as u32;
                let stop = pos as u32;
                if start < stop && (stop as usize) <= seq_len {
                    hits.push(TargetHit {
                        candidate: true,
                        seq: encoded_seq,
                        chrom_idx,
                        start,
                        stop,
                        strand: true,
                    });
                }
            }
        }
    }
    hits
}

/// Search 3prime orientation on reverse strand
pub fn search_3prime_reverse(
    chrom_idx: u32,
    seq: &[u8],
    pam_masks: &[u8],
    target_len: usize,
) -> Vec<TargetHit> {
    let mut hits = Vec::new();
    let pam_len = pam_masks.len();
    let flank_len = 2;
    let seq_len = seq.len();
    let total_window = target_len + pam_len + flank_len;
    if seq_len < total_window {
        return hits;
    }
    let rc_seq = dna::revcomp(seq);
    let min_rc_pos = target_len;
    let max_rc_pos = seq_len - pam_len - flank_len;
    for rc_pos in min_rc_pos..=max_rc_pos {
        if pam_matches_forward(pam_masks, &rc_seq, rc_pos) {
            let rc_target_start = rc_pos - target_len;
            let rc_window_end = rc_target_start + total_window;
            let window_bytes = &rc_seq[rc_target_start..rc_window_end];
            if let Ok(encoded_seq) = encode_2bit_u64(window_bytes) {
                let start_fwd = seq_len - rc_pos;
                let stop_fwd = seq_len - (rc_pos - target_len);
                if start_fwd < stop_fwd && stop_fwd <= seq_len {
                    hits.push(TargetHit {
                        candidate: true,
                        seq: encoded_seq,
                        chrom_idx,
                        start: start_fwd as u32,
                        stop: stop_fwd as u32,
                        strand: false,
                    });
                }
            }
        }
    }
    hits
}

/// Build Polars DataFrame from hits
pub fn build_dataframe(hits: &[TargetHit], chrom_names: &[String]) -> Result<DataFrame> {
    let mut candidate_vec = Vec::with_capacity(hits.len());
    let mut seq_vec = Vec::with_capacity(hits.len());
    let mut chrom_vec = Vec::with_capacity(hits.len());
    let mut start_vec = Vec::with_capacity(hits.len());
    let mut stop_vec = Vec::with_capacity(hits.len());
    let mut strand_vec = Vec::with_capacity(hits.len());

    for hit in hits {
        candidate_vec.push(hit.candidate);
        seq_vec.push(hit.seq);
        let name = chrom_names.get(hit.chrom_idx as usize).map(|s| s.as_str()).unwrap_or("unknown");
        chrom_vec.push(name);
        start_vec.push(hit.start);
        stop_vec.push(hit.stop);
        strand_vec.push(hit.strand);
    }

    let chrom_series = Series::new("chrom".into(), chrom_vec)
        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;

    let df = DataFrame::new(vec![
        Series::new("candidate".into(), candidate_vec).into(),
        Series::new("seq".into(), seq_vec).into(),
        chrom_series.into(),
        Series::new("start".into(), start_vec).into(),
        Series::new("stop".into(), stop_vec).into(),
        Series::new("strand".into(), strand_vec).into(),
    ])?;

    Ok(df)
}
