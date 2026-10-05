use hnsw_rs::prelude::*;

/// Calculate precomputed left-aligned target mask
#[inline(always)]
pub fn compute_target_mask(target_len: usize) -> u64 {
    let active_bits = 2 * target_len;
    if active_bits >= 64 {
        u64::MAX
    } else {
        !((1u64 << (64 - active_bits)) - 1)
    }
}

/// Calculate exact base Hamming distance between two 2-bit left-aligned u64 sequences using masked lane-collapse
#[inline(always)]
pub fn base_hamming_distance_masked(a: u64, b: u64, mask: u64) -> u32 {
    let x = a ^ b;
    (((x | (x >> 1)) & 0x5555_5555_5555_5555) & mask).count_ones()
}

/// Get segment base count and base offset for segment s out of m segments for target_len bases
#[inline(always)]
pub fn get_segment_bounds(s: usize, m: usize, target_len: usize) -> (usize, usize) {
    let base_len = target_len / m;
    let rem = target_len % m;

    let mut start_base = 0;
    for i in 0..s {
        start_base += base_len + if i < rem { 1 } else { 0 };
    }
    let seg_bases = base_len + if s < rem { 1 } else { 0 };
    (start_base, seg_bases)
}

/// Multi-Index Hashing (MIH) CSR Table for a single segment
#[derive(Debug)]
pub struct MihSegmentTable {
    pub shift: u32,
    pub mask: u64,
    pub offsets: Vec<u32>,
    pub target_ids: Vec<u32>,
}

impl MihSegmentTable {
    pub fn build(seq_vec: &[u64], target_len: usize, s: usize, m: usize) -> Self {
        let (start_base, seg_bases) = get_segment_bounds(s, m, target_len);
        let shift = (64 - 2 * (start_base + seg_bases)) as u32;
        let mask = if 2 * seg_bases >= 64 {
            u64::MAX
        } else {
            (1u64 << (2 * seg_bases)) - 1
        };

        let num_buckets = 1usize << (2 * seg_bases);
        let mut bucket_counts = vec![0u32; num_buckets];

        for &seq in seq_vec {
            let key = ((seq >> shift) & mask) as usize;
            bucket_counts[key] += 1;
        }

        let mut offsets = vec![0u32; num_buckets + 1];
        for i in 0..num_buckets {
            offsets[i + 1] = offsets[i] + bucket_counts[i];
        }

        let mut cursor = offsets.clone();
        let mut target_ids = vec![0u32; seq_vec.len()];

        for (i, &seq) in seq_vec.iter().enumerate() {
            let key = ((seq >> shift) & mask) as usize;
            let pos = cursor[key] as usize;
            target_ids[pos] = i as u32;
            cursor[key] += 1;
        }

        MihSegmentTable {
            shift,
            mask,
            offsets,
            target_ids,
        }
    }

    #[inline(always)]
    pub fn extract_key(&self, seq: u64) -> usize {
        ((seq >> self.shift) & self.mask) as usize
    }

    #[inline(always)]
    pub fn get_bucket_targets(&self, key: usize) -> &[u32] {
        if key + 1 < self.offsets.len() {
            let start = self.offsets[key] as usize;
            let end = self.offsets[key + 1] as usize;
            &self.target_ids[start..end]
        } else {
            &[]
        }
    }

    #[inline(always)]
    pub fn get_bucket_size(&self, key: usize) -> usize {
        if key + 1 < self.offsets.len() {
            (self.offsets[key + 1] - self.offsets[key]) as usize
        } else {
            0
        }
    }
}

/// DNA Hamming Distance metric for hnsw_rs
#[derive(Clone, Copy)]
pub struct DnaHammingDistance {
    pub mask: u64,
}

impl Distance<u64> for DnaHammingDistance {
    #[inline(always)]
    fn eval(&self, va: &[u64], vb: &[u64]) -> f32 {
        base_hamming_distance_masked(va[0], vb[0], self.mask) as f32
    }
}
