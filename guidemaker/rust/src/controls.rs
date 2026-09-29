use anyhow::Result;
use rand::RngExt; // Keeps the trait import for .gen() and .gen_bool()
use crate::base_hamming_distance_masked; // Pulls distance metric natively

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeedOrientation {
    ThreePrime, // e.g., Cas9 (Seed sits at the end of the guide before the PAM)
    FivePrime,  // e.g., Cas12a (Seed sits at the start of the guide after the PAM)
}


// Fast binary random sequence generation using target GC distribution bounds
#[inline(always)]
fn generate_random_2bit_seq(target_len: usize, gc_fraction: f64) -> u64 {
    let mut rng = rand::rng(); // Fetch the modern thread-local generator
    let mut val: u64 = 0;
    let gc_cutoff = (gc_fraction * u32::MAX as f64) as u32;

    for i in 0..target_len {
        let rand_val: u32 = rng.random(); 
        
        let code: u64 = if rand_val < gc_cutoff {
            if rng.random_bool(0.5) { 0b10 } else { 0b01 }
        } else {
            if rng.random_bool(0.5) { 0b00 } else { 0b11 }
        };
        
        let shift = 64 - 2 * (i + 1);
        val |= code << shift;
    }
    val
}




// Dynamically extract two 12-bit blocks (6 bases each) from the selected seed region
#[inline(always)]
fn extract_dynamic_seed_blocks(seq: u64, orientation: SeedOrientation, target_len: usize) -> [usize; 2] {
    match orientation {
        SeedOrientation::FivePrime => {
            let block0 = ((seq >> 52) & 0x0FFF) as usize;
            let block1 = ((seq >> 40) & 0x0FFF) as usize;
            [block0, block1]
        }
        SeedOrientation::ThreePrime => {
            let total_shift = 64 - (2 * target_len);
            let block0 = ((seq >> (total_shift + 12)) & 0x0FFF) as usize;
            let block1 = ((seq >> total_shift) & 0x0FFF) as usize;
            [block0, block1]
        }
    }
}

// Dynamically generate a bitmask covering exactly the seed bases
#[inline(always)]
fn compute_seed_mask(orientation: SeedOrientation, target_len: usize, seed_len: usize) -> u64 {
    let mut mask: u64 = 0;
    let base_mask = 0b11u64;

    match orientation {
        SeedOrientation::FivePrime => {
            for i in 0..seed_len {
                mask |= base_mask << (64 - 2 * (i + 1));
            }
        }
        SeedOrientation::ThreePrime => {
            let start_base = target_len - seed_len;
            for i in start_base..target_len {
                mask |= base_mask << (64 - 2 * (i + 1));
            }
        }
    }
    mask
}

pub fn generate_flexible_negative_controls(
    seq_vec: &[u64],
    orientation: SeedOrientation,
    target_len: usize,
    seed_len: usize,
    min_seed_mismatches: u32,
    gc_fraction: f64,
    requested_n: usize,
) -> Result<Vec<u64>> {
    let seed_mask = compute_seed_mask(orientation, target_len, seed_len);

    let mut seed_table0 = vec![Vec::<usize>::new(); 4096];
    let mut seed_table1 = vec![Vec::<usize>::new(); 4096];

    for (idx, &seq) in seq_vec.iter().enumerate() {
        let s_blocks = extract_dynamic_seed_blocks(seq, orientation, target_len);
        seed_table0[s_blocks[0]].push(idx);
        seed_table1[s_blocks[1]].push(idx);
    }

    let mut accepted_controls = Vec::with_capacity(requested_n);
    let mut iterations = 0;

    while accepted_controls.len() < requested_n && iterations < 1_000_000 {
        iterations += 1;
        let rand_seq = generate_random_2bit_seq(target_len, gc_fraction);
        let s_blocks = extract_dynamic_seed_blocks(rand_seq, orientation, target_len);

        if seed_table0[s_blocks[0]].is_empty() && seed_table1[s_blocks[1]].is_empty() {
            accepted_controls.push(rand_seq);
            continue;
        }

        let mut hit_toxic_overlap = false;
        for &ref_idx in seed_table0[s_blocks[0]].iter().chain(&seed_table1[s_blocks[1]]) {
            let ref_seq = seq_vec[ref_idx];
            let seed_dist = base_hamming_distance_masked(rand_seq, ref_seq, seed_mask);
            
            if seed_dist < min_seed_mismatches {
                hit_toxic_overlap = true;
                break;
            }
        }

        if !hit_toxic_overlap {
            accepted_controls.push(rand_seq);
        }
    }

    Ok(accepted_controls)
}
