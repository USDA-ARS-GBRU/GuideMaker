//! cfd.rs — Modular CFD scoring with YAML-loaded matrices + low-bit mismatch evaluator (LSR=20 nt)
//!
//! Exposes:
//!  - `CfdMatrix` struct: name + 20×(4×4) weights.
//!  - `load_cfd_matrix(path, name)` → CfdMatrix (validated).
//!  - `cfd_evaluate_lowbits(&CfdMatrix, tgt_low, cand_low)` → f32.
//!  - Helpers: mismatch bitmaps (low/high), precompute bounds, bounded evaluator.
//!
//! Bases (2-bit): A=0, C=1, G=2, T=3.

use anyhow::{Context, Result};
use serde::Deserialize;

pub const LSR_LEN: usize = 20;

// ---------- YAML schema & loader ----------

#[derive(Deserialize)]
struct CfdModelFile {
    models: Vec<CfdModelSpec>,
}

#[derive(Deserialize, Clone)]
struct CfdModelSpec {
    name: String,
    enzyme: Option<String>,
    version: Option<u32>,
    lsr_len: usize,                  // must be 20
    alphabet: [String; 4],           // ["A","C","G","T"]
    weights: Vec<[[f32; 4]; 4]>,     // length = 20, each is 4x4 tgt x cand
}

/// In-memory evaluator-ready matrix
#[derive(Clone)]
pub struct CfdMatrix {
    pub name: String,
    pub lsr_len: usize,                  // 20
    pub weights: Vec<[[f32; 4]; 4]>,     // 20 items
}

/// Load + validate + select by name
pub fn load_cfd_matrix(path: &std::path::Path, name: &str) -> Result<CfdMatrix> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read CFD YAML at {:?}", path))?;
    let parsed: CfdModelFile = serde_yaml::from_str(&text)
        .with_context(|| format!("Failed to parse CFD YAML at {:?}", path))?;

    let spec = parsed.models
        .iter()
        .find(|m| m.name == name)
        .ok_or_else(|| anyhow::anyhow!("CFD model '{}' not found in {:?}", name, path))?;

    // Basic validation
    if spec.lsr_len != 20 {
        anyhow::bail!("CFD model '{}' has lsr_len={}, expected 20", name, spec.lsr_len);
    }
    if spec.weights.len() != 20 {
        anyhow::bail!("CFD model '{}' must have 20 positions; found {}", name, spec.weights.len());
    }
    if spec.alphabet != ["A","C","G","T"] {
        anyhow::bail!("CFD model '{}' must use alphabet [A, C, G, T]", name);
    }

    Ok(CfdMatrix {
        name: spec.name.clone(),
        lsr_len: spec.lsr_len,
        weights: spec.weights.clone(),
    })
}

// ---------- mismatch-only evaluator (low-bit) ----------

const MASK_2BIT_LSB: u64 = 0x5555_5555_5555_5555; // ...0101 pattern: LSB of each 2-bit group

/// Mismatch bitmap for *low-bit* packing (positions in low 40 bits).
/// Each set bit indicates a mismatch at that position (bit i == position i).
#[inline(always)]
pub fn mismatch_bitmap_lowbits(wt: u64, off: u64, lsr_len_nt: usize) -> u64 {
    let mut x = wt ^ off;                 // any differing bit in a 2-bit symbol indicates mismatch
    x |= x >> 1;                          // fold 2-bit groups: set group LSB if either bit differs
    x & MASK_2BIT_LSB & ((1u64 << lsr_len_nt) - 1)
}

/// CFD score for *low-bit* packed sequences using the loaded matrix.
#[inline(always)]
pub fn cfd_evaluate_lowbits(
    cfd: &CfdMatrix,
    tgt_low: u64,
    cand_low: u64,
) -> f32 {
    // Build mismatch bitmap once
    let m = mismatch_bitmap_lowbits(tgt_low, cand_low, cfd.lsr_len);
    if m == 0 { return 1.0; }

    let mut score = 1.0f32;
    let mut bits = m as u32;
    while bits != 0 {
        let i = bits.trailing_zeros() as usize;         // position 0..19
        let t = ((tgt_low  >> (2*i)) & 0b11) as usize;  // 0..3
        let c = ((cand_low >> (2*i)) & 0b11) as usize;  // 0..3
        score *= cfd.weights[i][t][c];
        bits &= bits - 1;
    }
    score
}

// ---------- Optional helpers (for bounded evaluation / high-bit input) ----------

/// Mismatch bitmap for *high-bit* packing (positions occupy upper 40 bits).
#[inline(always)]
pub fn mismatch_bitmap_highbits(wt: u64, off: u64, lsr_len_nt: usize) -> u64 {
    let mut x = wt ^ off;
    x |= x >> 1;
    let top_mask = (!0u64) << (64 - 2 * lsr_len_nt); // keep only the upper 2*len bits
    x & MASK_2BIT_LSB & top_mask
}

/// Precompute per-position *maximum* mismatch weight across all base pairs.
/// Useful for bounded evaluation / early exit.
#[inline]
pub fn precompute_max_pos_weights_from_matrix(cfd: &CfdMatrix) -> [f32; LSR_LEN] {
    let mut out = [1.0f32; LSR_LEN];
    for pos in 0..LSR_LEN {
        let mut m = 0.0f32;
        for q in 0..4 {
            for t in 0..4 {
                if q == t { continue; }
                let w = cfd.weights[pos][q][t];
                if w > m { m = w; }
            }
        }
        out[pos] = if m > 0.0 { m } else { 1.0 };
    }
    out
}

/// Bounded CFD for *low-bit* packed sequences, using the loaded matrix.
#[inline(always)]
pub fn cfd_evaluate_lowbits_bounded(
    cfd: &CfdMatrix,
    tgt_low: u64,
    cand_low: u64,
    best_so_far: f32,
    max_weight_per_pos: &[f32; LSR_LEN],
) -> f32 {
    let m = mismatch_bitmap_lowbits(tgt_low, cand_low, cfd.lsr_len);

    // Upper bound of remaining multiplicative weights (lazy product over remaining positions).
    let mut remaining_prod_bound: f32 = 1.0;
    {
        let mut tmp = m as u32;
        while tmp != 0 {
            let pos = tmp.trailing_zeros() as usize;
            remaining_prod_bound *= max_weight_per_pos[pos];
            tmp &= tmp - 1;
        }
    }

    let mut score: f32 = 1.0;
    let mut mism = m as u32;

    while mism != 0 {
        if score * remaining_prod_bound <= best_so_far {
            break;
        }
        let pos = mism.trailing_zeros() as usize;

        // consume this position from bound (avoid div by 0; max is >= 1.0 typically)
        let this_max = max_weight_per_pos[pos];
        if this_max != 0.0 {
            remaining_prod_bound /= this_max;
        }

        let t = ((tgt_low  >> (2*pos)) & 0b11) as usize;
        let c = ((cand_low >> (2*pos)) & 0b11) as usize;
        score *= cfd.weights[pos][t][c];

        mism &= mism - 1;
    }

    score
}