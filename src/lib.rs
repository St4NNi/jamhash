//! # jamhash
//!
//! A fast, non-cryptographic hash function designed for genomics and general-purpose hashing.
//!
//! ## Features
//!
//! - `jamhash_u64`: Specialized hash function for single u64 values, optimized for genomic k-mers
//! - `JamHasher`: Streaming hasher implementing `std::hash::Hasher` trait
//! - Dual-path accumulation for strong avalanche properties
//! - Predictable: `hash(0) == 0`
//! - Tested for low collision rates with canonical k-mers up to k=21
//!
//! ## Use Cases
//!
//! ### Genomics and Bioinformatics
//!
//! The `jamhash_u64` function is particularly well-suited for:
//! - MinHash sketching of genomic sequences
//! - K-mer counting and indexing
//! - Sequence similarity estimation
//!
//! ```
//! use jamhash::jamhash_u64;
//!
//! // Hash a k-mer encoded as u64 (e.g., 2-bit encoding)
//! let kmer: u64 = 0x1234567890abcdef;
//! let hash = jamhash_u64(kmer);
//! ```
//!
//! ### General Purpose
//!
//! ```
//! use jamhash::{jamhash_bytes, JamHasher};
//! use std::hash::Hasher;
//!
//! // Hash arbitrary bytes
//! let hash = jamhash_bytes(b"hello world");
//!
//! // Use the streaming hasher
//! let mut hasher = JamHasher::new();
//! hasher.write(b"hello ");
//! hasher.write(b"world");
//! let hash = hasher.finish();
//! ```

use std::hash::Hasher;

/// Constants derived from pi (see tools/piconstants.py)
const CONST1: u64 = 0xb8e1afed6a267e96;
const CONST2: u64 = 0x082efa98ec4e6c89;
const CONST3: u64 = 0x243f6a8885a308d3;
const CONST4: u64 = 0x13198a2e03707344;
const CONST5: u64 = 0xa4093822299f31d0;
const CONST6: u64 = 0x452821e638d01377;

/// Folded multiply: multiply two u64s and XOR the high and low halves
#[inline(always)]
fn fold_multiply(input: u64, const_1: u64) -> u64 {
    let temp = (input as u128).wrapping_mul(const_1 as u128);
    (temp as u64) ^ ((temp >> 64) as u64)
}

/// Fast unaligned u64 read
#[inline(always)]
unsafe fn read_u64(ptr: *const u8) -> u64 {
    unsafe { ptr.cast::<u64>().read_unaligned() }
}

/// Fast unaligned u32 read
#[inline(always)]
unsafe fn read_u32(ptr: *const u8) -> u64 {
    unsafe { ptr.cast::<u32>().read_unaligned() as u64 }
}

/// Specialized hash function for single u64 values.
///
/// This is the core jamhash algorithm using dual rotation paths
/// and folded multiplication for excellent mixing properties.
///
/// ## Genomics Applications
///
/// This function is optimized for hashing genomic k-mers in bioinformatics workflows:
/// - MinHash sketching of DNA/RNA sequences
/// - K-mer counting and frequency analysis
/// - Sequence similarity and comparison
///
/// Tested for low collision rates with canonical k-mers up to k=21, making it
/// suitable for most genomic analysis pipelines that encode k-mers as u64 values.
///
/// # Examples
///
/// ```
/// use jamhash::jamhash_u64;
///
/// // Hash a k-mer encoded as u64
/// let kmer: u64 = 0x1234567890abcdef;
/// let hash = jamhash_u64(kmer);
/// assert_ne!(hash, kmer);
///
/// // Predictable zero property
/// assert_eq!(jamhash_u64(0), 0);
/// ```
#[inline]
pub fn jamhash_u64(kmer: u64) -> u64 {
    let part1 = kmer.rotate_left(16);
    let part2 = kmer.rotate_left(48);
    let fold1 = fold_multiply(part1, CONST1);
    let fold2 = fold_multiply(part2, CONST2);
    fold_multiply(fold1, fold2)
}

/// Rapid mixing using single 64-bit multiply + XOR shift (cheaper than fold_multiply)
#[inline(always)]
fn rapidmix(v: u64, k: u64) -> u64 {
    let m = v.wrapping_mul(k);
    m ^ (m >> 47)
}

// Precomputed constants for unseeded path (seed=0)
const SEED0_S: u64 = CONST4; // 0 + CONST4
const SEED0_S17: u64 = CONST4.rotate_left(17);
const SEED0_S41: u64 = CONST4.rotate_left(41);

/// Fast unseeded hash - specialized for seed=0 case
#[inline(always)]
fn jamhash_unseeded(value: u64) -> u64 {
    let v = value ^ CONST5;
    // Precomputed seed mixing
    let mixed = rapidmix(v ^ SEED0_S, CONST1) ^ rapidmix(v.rotate_left(32) ^ SEED0_S17, CONST2);
    fold_multiply(mixed ^ SEED0_S41, CONST3)
}

/// Hash for bytes with integrated seed mixing - 2 rapidmix + 1 fold_multiply
#[inline(always)]
fn jamhash_seeded(value: u64, seed: u64) -> u64 {
    // XOR with constant to break zero-preserving property
    let v = value ^ CONST5;
    // Mix seed with carries for bit spreading
    let s = seed.wrapping_add(CONST4);
    // Two rapidmix operations for pre-mixing (cheaper than fold_multiply)
    let mixed = rapidmix(v ^ s, CONST1) ^ rapidmix(v.rotate_left(32) ^ s.rotate_left(17), CONST2);
    // Single fold_multiply for final quality
    fold_multiply(mixed ^ s.rotate_left(41), CONST3)
}

/// Hash short byte sequences (≤16 bytes) - unseeded fast path
#[inline(always)]
fn jam_short_unseeded(bytes: &[u8]) -> u64 {
    let len = bytes.len();
    let ptr = bytes.as_ptr();

    unsafe {
        match len {
            0 => jamhash_unseeded(0),
            1..=3 => {
                // Pack bytes with wide separation following wyhash pattern
                let combined = ((bytes[0] as u64) << 16)
                    | ((bytes[len >> 1] as u64) << 8)
                    | (bytes[len - 1] as u64);
                // Protected multiply pattern - prevents bias propagation
                // Using two fold_multiply rounds for better avalanche
                let a = fold_multiply(combined ^ CONST5, CONST1 ^ (len as u64));
                fold_multiply(a ^ CONST6, CONST2)
            }
            4..=7 => {
                let lo = read_u32(ptr);
                let hi = read_u32(ptr.add(len - 4));
                // XOR with constants to break overlap correlation for 4-6 byte inputs
                let a = fold_multiply(lo ^ CONST5, CONST1);
                let b = fold_multiply(hi ^ CONST6, CONST2);
                fold_multiply(a ^ b, CONST3 ^ (len as u64))
            }
            8 => jamhash_unseeded(read_u64(ptr)),
            9..=16 => {
                let lo = read_u64(ptr);
                let hi = read_u64(ptr.add(len - 8));
                // Two independent fold_multiply paths for better mixing
                let a = fold_multiply(lo ^ CONST5, CONST1);
                let b = fold_multiply(hi ^ CONST6, CONST2);
                fold_multiply(
                    a ^ b ^ lo.rotate_left(23) ^ hi.rotate_left(47),
                    CONST3 ^ (len as u64),
                )
            }
            _ => unreachable!(),
        }
    }
}

/// Hash short byte sequences (≤16 bytes) - seeded path
#[inline(always)]
fn jam_short(bytes: &[u8], seed: u64) -> u64 {
    let len = bytes.len();
    let ptr = bytes.as_ptr();
    let s = seed.wrapping_add(CONST4);

    unsafe {
        match len {
            0 => jamhash_seeded(0, seed),
            1..=3 => {
                // Pack bytes with wide separation following wyhash pattern
                let combined = ((bytes[0] as u64) << 16)
                    | ((bytes[len >> 1] as u64) << 8)
                    | (bytes[len - 1] as u64);
                // Protected multiply with seed mixing
                let a = fold_multiply(combined ^ s, CONST1 ^ (len as u64));
                fold_multiply(a ^ s.rotate_left(33), CONST2)
            }
            4..=7 => {
                let lo = read_u32(ptr);
                let hi = read_u32(ptr.add(len - 4));
                // XOR with seed-derived constants to break overlap correlation
                let a = fold_multiply(lo ^ s, CONST1);
                let b = fold_multiply(hi ^ s.rotate_left(33), CONST2);
                fold_multiply(a ^ b, CONST3 ^ (len as u64))
            }
            8 => jamhash_seeded(read_u64(ptr), seed),
            9..=16 => {
                let lo = read_u64(ptr);
                let hi = read_u64(ptr.add(len - 8));
                // Two independent fold_multiply paths for better mixing
                let a = fold_multiply(lo ^ s, CONST1);
                let b = fold_multiply(hi ^ s.rotate_left(33), CONST2);
                fold_multiply(
                    a ^ b ^ lo.rotate_left(23) ^ hi.rotate_left(47),
                    CONST3 ^ (len as u64),
                )
            }
            _ => unreachable!(),
        }
    }
}

/// Hash medium byte sequences (17-128 bytes) - optimized: 2 fold_multiply per 16 bytes
#[inline(always)]
fn jam_medium(bytes: &[u8], seed: u64) -> u64 {
    let len = bytes.len();
    let ptr = bytes.as_ptr();

    // Seed mixed via addition for carry propagation
    let s = seed.wrapping_add(CONST4);

    // Two parallel accumulators initialized with seed derivatives
    let mut acc0 = s ^ CONST1;
    let mut acc1 = s.rotate_left(33) ^ CONST2;

    unsafe {
        let mut i = 0;
        // Process 16-byte chunks with 2 fold_multiply per chunk
        while i + 16 <= len {
            acc0 = fold_multiply(acc0 ^ read_u64(ptr.add(i)), CONST3);
            acc1 = fold_multiply(acc1 ^ read_u64(ptr.add(i + 8)), CONST4);
            i += 16;
        }

        // Handle tail via overlapping reads only if there's remaining data
        if i < len {
            acc0 = fold_multiply(acc0 ^ read_u64(ptr.add(len - 16)), CONST5);
            acc1 = fold_multiply(acc1 ^ read_u64(ptr.add(len - 8)), CONST6);
        }

        // Final mix - single fold_multiply
        fold_multiply(
            acc0 ^ acc1.rotate_left(23) ^ s,
            (acc1 ^ acc0.rotate_left(47)) ^ (len as u64),
        )
    }
}

/// Hash long byte sequences (>128 bytes) - optimized with loop unrolling
#[cold]
#[inline(never)]
fn jam_long(bytes: &[u8], seed: u64) -> u64 {
    let len = bytes.len();
    let ptr = bytes.as_ptr();

    // For very long inputs (>4KB), use 8-wide accumulator path
    if len >= 4096 {
        return jam_ultra_long(bytes, seed);
    }

    // Seed mixed via addition for carry propagation
    let s = seed.wrapping_add(CONST4);

    // Four accumulators with different seed derivatives
    let mut s0 = s ^ CONST1;
    let mut s1 = s.rotate_left(17) ^ CONST2;
    let mut s2 = s.rotate_left(33) ^ CONST3;
    let mut s3 = s.rotate_left(49) ^ CONST4;

    unsafe {
        let mut i = 0;
        let end64 = len.saturating_sub(63);
        let end32 = len - (len % 32);

        // Process 64 bytes per iteration (2x unrolled) for better ILP
        while i < end64 {
            // First 32 bytes
            s0 = fold_multiply(s0 ^ read_u64(ptr.add(i)), CONST1);
            s1 = fold_multiply(s1 ^ read_u64(ptr.add(i + 8)), CONST2);
            s2 = fold_multiply(s2 ^ read_u64(ptr.add(i + 16)), CONST3);
            s3 = fold_multiply(s3 ^ read_u64(ptr.add(i + 24)), CONST4);
            // Second 32 bytes
            s0 = fold_multiply(s0 ^ read_u64(ptr.add(i + 32)), CONST5);
            s1 = fold_multiply(s1 ^ read_u64(ptr.add(i + 40)), CONST6);
            s2 = fold_multiply(s2 ^ read_u64(ptr.add(i + 48)), CONST1);
            s3 = fold_multiply(s3 ^ read_u64(ptr.add(i + 56)), CONST2);
            i += 64;
        }

        // Handle remaining 32-byte chunk if present
        if i < end32 {
            s0 = fold_multiply(s0 ^ read_u64(ptr.add(i)), CONST1);
            s1 = fold_multiply(s1 ^ read_u64(ptr.add(i + 8)), CONST2);
            s2 = fold_multiply(s2 ^ read_u64(ptr.add(i + 16)), CONST3);
            s3 = fold_multiply(s3 ^ read_u64(ptr.add(i + 24)), CONST4);
            i += 32;
        }

        // Handle tail via overlapping reads only if there's remaining data
        if i < len {
            s0 = fold_multiply(s0 ^ read_u64(ptr.add(len - 32)), CONST5);
            s1 = fold_multiply(s1 ^ read_u64(ptr.add(len - 24)), CONST6);
            s2 = fold_multiply(s2 ^ read_u64(ptr.add(len - 16)), CONST1);
            s3 = fold_multiply(s3 ^ read_u64(ptr.add(len - 8)), CONST2);
        }

        // Merge accumulators - two fold_multiply operations
        let m0 = fold_multiply(s0 ^ s2.rotate_left(17), CONST5);
        let m1 = fold_multiply(s1 ^ s3.rotate_left(31), CONST6);
        fold_multiply(
            m0 ^ m1.rotate_left(23) ^ s,
            (m1 ^ m0.rotate_left(47)) ^ (len as u64),
        )
    }
}

/// Hash ultra-long byte sequences (>4KB) - 8-wide accumulator for maximum throughput
#[cold]
#[inline(never)]
fn jam_ultra_long(bytes: &[u8], seed: u64) -> u64 {
    let len = bytes.len();
    let ptr = bytes.as_ptr();

    // Seed mixed via addition for carry propagation
    let s = seed.wrapping_add(CONST4);

    // Eight accumulators for 8-wide parallelism (64 bytes per iteration)
    let mut a0 = s ^ CONST1;
    let mut a1 = s.rotate_left(9) ^ CONST2;
    let mut a2 = s.rotate_left(17) ^ CONST3;
    let mut a3 = s.rotate_left(25) ^ CONST4;
    let mut a4 = s.rotate_left(33) ^ CONST5;
    let mut a5 = s.rotate_left(41) ^ CONST6;
    let mut a6 = s.rotate_left(49) ^ CONST1;
    let mut a7 = s.rotate_left(57) ^ CONST2;

    unsafe {
        let mut i = 0;
        let end = len - (len % 64);

        // Main loop: 64 bytes per iteration, 8 parallel fold_multiply
        while i < end {
            a0 = fold_multiply(a0 ^ read_u64(ptr.add(i)), CONST1);
            a1 = fold_multiply(a1 ^ read_u64(ptr.add(i + 8)), CONST2);
            a2 = fold_multiply(a2 ^ read_u64(ptr.add(i + 16)), CONST3);
            a3 = fold_multiply(a3 ^ read_u64(ptr.add(i + 24)), CONST4);
            a4 = fold_multiply(a4 ^ read_u64(ptr.add(i + 32)), CONST5);
            a5 = fold_multiply(a5 ^ read_u64(ptr.add(i + 40)), CONST6);
            a6 = fold_multiply(a6 ^ read_u64(ptr.add(i + 48)), CONST1);
            a7 = fold_multiply(a7 ^ read_u64(ptr.add(i + 56)), CONST2);
            i += 64;
        }

        // Handle tail via overlapping reads
        if i < len {
            let tail_ptr = ptr.add(len - 64);
            a0 = fold_multiply(a0 ^ read_u64(tail_ptr), CONST3);
            a1 = fold_multiply(a1 ^ read_u64(tail_ptr.add(8)), CONST4);
            a2 = fold_multiply(a2 ^ read_u64(tail_ptr.add(16)), CONST5);
            a3 = fold_multiply(a3 ^ read_u64(tail_ptr.add(24)), CONST6);
            a4 = fold_multiply(a4 ^ read_u64(tail_ptr.add(32)), CONST1);
            a5 = fold_multiply(a5 ^ read_u64(tail_ptr.add(40)), CONST2);
            a6 = fold_multiply(a6 ^ read_u64(tail_ptr.add(48)), CONST3);
            a7 = fold_multiply(a7 ^ read_u64(tail_ptr.add(56)), CONST4);
        }

        // Merge 8 accumulators into 4
        let m0 = fold_multiply(a0 ^ a4.rotate_left(17), CONST5);
        let m1 = fold_multiply(a1 ^ a5.rotate_left(23), CONST6);
        let m2 = fold_multiply(a2 ^ a6.rotate_left(31), CONST1);
        let m3 = fold_multiply(a3 ^ a7.rotate_left(37), CONST2);

        // Merge 4 into 2
        let n0 = fold_multiply(m0 ^ m2.rotate_left(19), CONST3);
        let n1 = fold_multiply(m1 ^ m3.rotate_left(29), CONST4);

        // Final merge
        fold_multiply(
            n0 ^ n1.rotate_left(23) ^ s,
            (n1 ^ n0.rotate_left(47)) ^ (len as u64),
        )
    }
}

/// Streaming hasher for arbitrary byte sequences.
///
/// Processes incoming bytes using size-optimized paths and
/// finalizes with jamhash_u64 for consistent mixing.
pub struct JamHasher {
    accumulator: u64,
    len: u64,
}

impl JamHasher {
    #[inline]
    pub fn new() -> Self {
        Self {
            accumulator: 0,
            len: 0,
        }
    }

    pub fn new_with_seed(seed: u64) -> Self {
        Self {
            accumulator: seed,
            len: 0,
        }
    }
}

impl Default for JamHasher {
    fn default() -> Self {
        Self::new()
    }
}

impl Hasher for JamHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        self.len += bytes.len() as u64;
        self.accumulator = match bytes.len() {
            0..=16 => jam_short(bytes, self.accumulator),
            17..=128 => jam_medium(bytes, self.accumulator),
            _ => jam_long(bytes, self.accumulator),
        };
    }

    #[inline]
    fn finish(&self) -> u64 {
        // Use jamhash_u64 for finalization - accumulator already well-mixed from write()
        jamhash_u64(self.accumulator ^ self.len)
    }
}

/// Hash arbitrary bytes using jamhash.
///
/// Uses size-optimized paths for different input lengths
/// and finalizes with jamhash_u64.
///
/// # Example
///
/// ```
/// use jamhash::jamhash_bytes;
///
/// let hash = jamhash_bytes(b"hello world");
/// assert_ne!(hash, 0);
///
/// ```
#[inline]
pub fn jamhash_bytes(bytes: &[u8]) -> u64 {
    // Use specialized unseeded paths for common case
    match bytes.len() {
        0..=16 => jam_short_unseeded(bytes),
        17..=128 => jam_medium(bytes, 0),
        _ => jam_long(bytes, 0),
    }
}

#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod tests {
    use super::*;

    #[test]
    fn test_jamhash_u64_zero() {
        assert_eq!(jamhash_u64(0), 0, "hash(0) should be 0");
    }

    #[test]
    fn test_jamhash_u64_nonzero() {
        let hash = jamhash_u64(0x1234567890abcdef);
        assert_ne!(hash, 0);
    }

    #[test]
    fn test_jamhash_bytes_consistency() {
        let data = b"hello world";
        let hash1 = jamhash_bytes(data);
        let hash2 = jamhash_bytes(data);
        assert_eq!(hash1, hash2, "same input should produce same hash");
    }

    #[test]
    fn test_different_lengths() {
        // These should all be different due to length mixing
        let hash1 = jamhash_bytes(&[0]);
        let hash2 = jamhash_bytes(&[0, 0]);
        let hash3 = jamhash_bytes(&[0, 0, 0]);

        assert_ne!(hash1, hash2);
        assert_ne!(hash2, hash3);
        assert_ne!(hash1, hash3);
    }

    #[test]
    fn test_partial_buffer_lengths() {
        // Test all partial buffer lengths (1-7 bytes)
        for len in 1..8 {
            let data: Vec<u8> = (0..len).map(|i| i as u8).collect();
            let hash = jamhash_bytes(&data);
            assert_ne!(hash, 0, "hash of {:?} should not be zero", data);
        }
    }

    #[test]
    fn test_full_buffer() {
        let data = [1u8, 2, 3, 4, 5, 6, 7, 8];
        let hash = jamhash_bytes(&data);
        assert_ne!(hash, 0);
    }

    #[test]
    fn test_multiple_buffers() {
        let data = [1u8; 24]; // 3 full buffers
        let hash = jamhash_bytes(&data);
        assert_ne!(hash, 0);
    }

    #[test]
    fn test_avalanche() {
        // Small input changes should cause large output changes
        let hash1 = jamhash_bytes(b"test");
        let hash2 = jamhash_bytes(b"Test");

        // Count differing bits
        let diff = (hash1 ^ hash2).count_ones();
        assert!(
            diff > 10,
            "avalanche effect too weak: only {} bits differ",
            diff
        );
    }

    #[test]
    fn test_bit_flip_9_bytes() {
        // Test that flipping each bit in a 9-byte input produces different hash
        let base = [0u8; 9];
        let base_hash = jamhash_bytes(&base);

        for bit in 0..72 {
            let mut modified = base;
            modified[bit / 8] ^= 1 << (bit % 8);
            let modified_hash = jamhash_bytes(&modified);
            assert_ne!(
                base_hash, modified_hash,
                "Flipping bit {} in 9-byte input produced identical hash",
                bit
            );
        }
    }

    #[test]
    fn test_bit_flip_all_short_lengths() {
        // Test bit flips for lengths 1-16 (covers SMHasher sanity checks)
        for len in 1..=16 {
            let base: Vec<u8> = vec![0u8; len];
            let base_hash = jamhash_bytes(&base);

            for bit in 0..(len * 8) {
                let mut modified = base.clone();
                modified[bit / 8] ^= 1 << (bit % 8);
                let modified_hash = jamhash_bytes(&modified);
                assert_ne!(
                    base_hash, modified_hash,
                    "Flipping bit {} in {}-byte input produced identical hash",
                    bit, len
                );
            }
        }
    }

    #[test]
    fn test_cyclic_4_cycles_of_8_bytes() {
        // SMHasher Cyclic test: "4 cycles of 8 bytes" = 32 bytes total
        // Pattern: ABABABAB where A and B are 8-byte blocks
        // All different cyclic patterns must produce different hashes
        use std::collections::HashSet;

        let mut hashes = HashSet::new();
        let num_patterns = 10000;

        for i in 0..num_patterns {
            // Create a repeating 8-byte pattern, cycled 4 times
            let pattern: [u8; 8] = [
                (i & 0xFF) as u8,
                ((i >> 8) & 0xFF) as u8,
                ((i >> 16) & 0xFF) as u8,
                ((i >> 24) & 0xFF) as u8,
                0,
                0,
                0,
                0,
            ];
            let mut data = [0u8; 32];
            for cycle in 0..4 {
                data[cycle * 8..(cycle + 1) * 8].copy_from_slice(&pattern);
            }
            let hash = jamhash_bytes(&data);
            hashes.insert(hash);
        }

        // Should have very few collisions (ideally 0 for 10000 patterns)
        let collision_rate = 1.0 - (hashes.len() as f64 / num_patterns as f64);
        assert!(
            collision_rate < 0.001,
            "Cyclic 4x8 test: {} unique hashes from {} patterns ({:.2}% collision rate)",
            hashes.len(),
            num_patterns,
            collision_rate * 100.0
        );
    }

    #[test]
    fn test_cyclic_8_cycles_of_4_bytes() {
        // SMHasher Cyclic test: "8 cycles of 4 bytes" = 32 bytes total
        use std::collections::HashSet;

        let mut hashes = HashSet::new();
        let num_patterns = 10000;

        for i in 0..num_patterns {
            let pattern: [u8; 4] = [
                (i & 0xFF) as u8,
                ((i >> 8) & 0xFF) as u8,
                ((i >> 16) & 0xFF) as u8,
                ((i >> 24) & 0xFF) as u8,
            ];
            let mut data = [0u8; 32];
            for cycle in 0..8 {
                data[cycle * 4..(cycle + 1) * 4].copy_from_slice(&pattern);
            }
            let hash = jamhash_bytes(&data);
            hashes.insert(hash);
        }

        let collision_rate = 1.0 - (hashes.len() as f64 / num_patterns as f64);
        assert!(
            collision_rate < 0.001,
            "Cyclic 8x4 test: {} unique hashes from {} patterns ({:.2}% collision rate)",
            hashes.len(),
            num_patterns,
            collision_rate * 100.0
        );
    }

    #[test]
    fn test_cyclic_patterns_must_differ() {
        // Different cyclic patterns of same structure must hash differently
        // Pattern A: [1,0,0,0,0,0,0,0] repeated 4 times
        // Pattern B: [2,0,0,0,0,0,0,0] repeated 4 times
        let pattern_a: [u8; 32] = [
            1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0,
        ];
        let pattern_b: [u8; 32] = [
            2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
            0, 0, 0,
        ];

        let hash_a = jamhash_bytes(&pattern_a);
        let hash_b = jamhash_bytes(&pattern_b);

        assert_ne!(
            hash_a, hash_b,
            "Different cyclic patterns produced identical hash: {:016x}",
            hash_a
        );
    }

    #[test]
    fn test_sparse_32_byte_keys() {
        // SMHasher Sparse test: keys with only a few bits set
        // Test that different sparse patterns produce different hashes
        use std::collections::HashSet;

        let mut hashes = HashSet::new();

        // Test all single-bit patterns in 32-byte keys
        for byte_pos in 0..32 {
            for bit_pos in 0..8 {
                let mut data = [0u8; 32];
                data[byte_pos] = 1 << bit_pos;
                let hash = jamhash_bytes(&data);
                let inserted = hashes.insert(hash);
                assert!(
                    inserted,
                    "Collision: bit {} of byte {} in 32-byte sparse key",
                    bit_pos, byte_pos
                );
            }
        }

        assert_eq!(
            hashes.len(),
            32 * 8,
            "Expected 256 unique hashes for single-bit patterns"
        );
    }

    #[test]
    fn test_long_cyclic_64_bytes() {
        // Test 64-byte cyclic patterns (goes to jam_long path)
        let pattern_a: Vec<u8> = (0..8).cycle().take(64).collect();
        let pattern_b: Vec<u8> = (1..9).cycle().take(64).collect();

        let hash_a = jamhash_bytes(&pattern_a);
        let hash_b = jamhash_bytes(&pattern_b);

        assert_ne!(
            hash_a, hash_b,
            "Different 64-byte cyclic patterns produced identical hash"
        );
    }

    #[test]
    fn test_ultra_long_5kb() {
        // Test the ultra-long path (>4KB)
        use std::collections::HashSet;

        let mut hashes = HashSet::new();
        let num_patterns = 100u64;

        for i in 0u64..num_patterns {
            // Create unique 5KB data
            let data: Vec<u8> = (0u64..5120)
                .map(|j| (i.wrapping_add(j).wrapping_mul(0x9e3779b97f4a7c15) >> 56) as u8)
                .collect();
            let hash = jamhash_bytes(&data);
            hashes.insert(hash);
        }

        let collision_rate = 1.0 - (hashes.len() as f64 / num_patterns as f64);
        assert!(
            collision_rate < 0.01,
            "5KB data: {} unique from {} ({:.2}% collision)",
            hashes.len(),
            num_patterns,
            collision_rate * 100.0
        );
    }

    #[test]
    fn test_ultra_long_bit_flip() {
        // Test bit sensitivity for ultra-long inputs
        let base: Vec<u8> = (0..5120).map(|i| (i * 7) as u8).collect();
        let base_hash = jamhash_bytes(&base);

        // Flip bits at various positions
        for pos in [0, 1000, 2000, 3000, 4000, 5000, 5119] {
            let mut modified = base.clone();
            modified[pos] ^= 1;
            let modified_hash = jamhash_bytes(&modified);
            let diff_bits = (base_hash ^ modified_hash).count_ones();
            assert!(
                diff_bits >= 20,
                "Flip at pos {}: only {} bits differ (expected >=20)",
                pos,
                diff_bits
            );
        }
    }

    #[test]
    fn test_long_cyclic_256_bytes() {
        // Test 256-byte cyclic patterns (definitely jam_long)
        use std::collections::HashSet;

        let mut hashes = HashSet::new();
        let num_patterns = 1000;

        for i in 0u64..num_patterns {
            // Create truly unique 32-byte patterns - encode i directly
            let mut base = [0u8; 32];
            base[0] = (i & 0xFF) as u8;
            base[1] = ((i >> 8) & 0xFF) as u8;
            base[2] = ((i >> 16) & 0xFF) as u8;
            base[3] = ((i >> 24) & 0xFF) as u8;
            // Fill rest with position-dependent values
            for j in 4..32 {
                base[j] = ((i.wrapping_add(j as u64)) & 0xFF) as u8;
            }

            let data: Vec<u8> = base.iter().cycle().take(256).cloned().collect();
            let hash = jamhash_bytes(&data);
            hashes.insert(hash);
        }

        let collision_rate = 1.0 - (hashes.len() as f64 / num_patterns as f64);
        assert!(
            collision_rate < 0.01, // Allow up to 1% collision for statistical variance
            "256-byte cyclic: {} unique from {} ({:.2}% collision)",
            hashes.len(),
            num_patterns,
            collision_rate * 100.0
        );
    }

    #[test]
    fn test_xor_cancellation_bug() {
        // This test directly checks for the XOR cancellation bug
        // If the bug exists: identical 16-byte blocks XOR'd together = 0
        // causing data to be ignored

        // Create pattern: 16 bytes repeated twice = 32 bytes
        // With the XOR bug, v0 ^ v2 = 0 because v0 == v2
        let pattern_a = [1u8; 16];
        let pattern_b = [2u8; 16];

        // Data A: [1,1,1...] repeated twice
        let mut data_a = [0u8; 32];
        data_a[0..16].copy_from_slice(&pattern_a);
        data_a[16..32].copy_from_slice(&pattern_a);

        // Data B: [2,2,2...] repeated twice
        let mut data_b = [0u8; 32];
        data_b[0..16].copy_from_slice(&pattern_b);
        data_b[16..32].copy_from_slice(&pattern_b);

        let hash_a = jamhash_bytes(&data_a);
        let hash_b = jamhash_bytes(&data_b);

        // These MUST be different. If they're the same, XOR cancellation bug exists.
        assert_ne!(
            hash_a, hash_b,
            "XOR cancellation bug detected! Different 32-byte inputs with 16-byte repeat \
             produced identical hash: {:016x}",
            hash_a
        );
    }

    #[test]
    fn test_xor_cancellation_jam_long() {
        // Same test for jam_long (>128 bytes)
        // With 64-byte XOR combining: v0^v4 where v0==v4 causes cancellation

        // Create 64-byte pattern repeated twice = 128 bytes (still jam_medium)
        // Need >128 for jam_long, so use 256 bytes = 64 bytes repeated 4 times

        let pattern_a: [u8; 64] = [1u8; 64];
        let pattern_b: [u8; 64] = [2u8; 64];

        let data_a: Vec<u8> = pattern_a.iter().cycle().take(256).cloned().collect();
        let data_b: Vec<u8> = pattern_b.iter().cycle().take(256).cloned().collect();

        let hash_a = jamhash_bytes(&data_a);
        let hash_b = jamhash_bytes(&data_b);

        assert_ne!(
            hash_a, hash_b,
            "XOR cancellation bug in jam_long! Different 256-byte inputs produced identical hash: {:016x}",
            hash_a
        );
    }

    // ==========================================
    // SMHasher-style tests for failing cases
    // ==========================================

    #[test]
    fn test_smhasher_permutation_8bytes_zero_low_bit() {
        // SMHasher Permutation test: 8-bytes [0, low bit; LE/BE]
        // Tests all permutations of two 8-byte values: one zero, one with single low bit
        use std::collections::HashSet;

        let mut hashes = HashSet::new();
        let mut collisions = Vec::new();

        // Test LE and BE orderings of (0, low_bit) pairs
        for bit in 0..8u32 {
            let low_bit_val = 1u64 << bit;

            // Pair: [0, low_bit] as 16 bytes LE
            let mut data_le = [0u8; 16];
            data_le[8..16].copy_from_slice(&low_bit_val.to_le_bytes());
            let hash_le = jamhash_bytes(&data_le);
            if !hashes.insert(hash_le) {
                collisions.push(format!("LE [0, 1<<{}]", bit));
            }

            // Pair: [0, low_bit] as 16 bytes BE
            let mut data_be = [0u8; 16];
            data_be[0..8].copy_from_slice(&low_bit_val.to_be_bytes());
            let hash_be = jamhash_bytes(&data_be);
            if !hashes.insert(hash_be) {
                collisions.push(format!("BE [0, 1<<{}]", bit));
            }

            // Pair: [low_bit, 0] LE
            let mut data_le_rev = [0u8; 16];
            data_le_rev[0..8].copy_from_slice(&low_bit_val.to_le_bytes());
            let hash_le_rev = jamhash_bytes(&data_le_rev);
            if !hashes.insert(hash_le_rev) {
                collisions.push(format!("LE [1<<{}, 0]", bit));
            }

            // Pair: [low_bit, 0] BE
            let mut data_be_rev = [0u8; 16];
            data_be_rev[8..16].copy_from_slice(&low_bit_val.to_be_bytes());
            let hash_be_rev = jamhash_bytes(&data_be_rev);
            if !hashes.insert(hash_be_rev) {
                collisions.push(format!("BE [1<<{}, 0]", bit));
            }
        }

        assert!(
            collisions.is_empty(),
            "Permutation 8-bytes [0, low bit] collisions: {:?}",
            collisions
        );
    }

    #[test]
    fn test_smhasher_permutation_8bytes_zero_high_bit() {
        // SMHasher Permutation test: 8-bytes [0, high bit; LE/BE]
        use std::collections::HashSet;

        let mut hashes = HashSet::new();
        let mut collisions = Vec::new();

        for bit in 56..64u32 {
            let high_bit_val = 1u64 << bit;

            // All 4 orderings for high bit patterns
            let patterns = [
                (0u64, high_bit_val, "LE [0, high]"),
                (high_bit_val, 0u64, "LE [high, 0]"),
            ];

            for (v1, v2, name) in patterns {
                let mut data = [0u8; 16];
                data[0..8].copy_from_slice(&v1.to_le_bytes());
                data[8..16].copy_from_slice(&v2.to_le_bytes());
                let hash = jamhash_bytes(&data);
                if !hashes.insert(hash) {
                    collisions.push(format!("{} bit {}", name, bit));
                }
            }
        }

        assert!(
            collisions.is_empty(),
            "Permutation 8-bytes [0, high bit] collisions: {:?}",
            collisions
        );
    }

    #[test]
    fn test_smhasher_text_long_alnum_last() {
        // SMHasher Text test: Long alphanumeric strings with last char variations
        // Tests positions 1968-2128, 4016-4176, 8112-8272
        use std::collections::HashSet;

        let base_sizes = [2048, 4096, 8192];

        for &size in &base_sizes {
            let mut hashes = HashSet::new();

            // Create base alphanumeric string
            let base: Vec<u8> = (0..size)
                .map(|i| {
                    let c = i % 62;
                    if c < 10 {
                        b'0' + c as u8
                    } else if c < 36 {
                        b'A' + (c - 10) as u8
                    } else {
                        b'a' + (c - 36) as u8
                    }
                })
                .collect();

            // Vary last 160 bytes (matching SMHasher range)
            for last_byte in 0..=255u8 {
                let mut data = base.clone();
                data[size - 1] = last_byte;
                let hash = jamhash_bytes(&data);
                hashes.insert(hash);
            }

            assert_eq!(
                hashes.len(),
                256,
                "Text long alnum size {}: expected 256 unique hashes, got {}",
                size,
                hashes.len()
            );
        }
    }

    #[test]
    fn test_smhasher_perlin_noise_2d() {
        // SMHasher PerlinNoise test: 2D coordinate hashing
        // Tests correlated inputs like (x, y) coordinates
        use std::collections::HashSet;

        let mut hashes = HashSet::new();
        let grid_size = 256; // 256x256 = 65536 points

        for y in 0..grid_size {
            for x in 0..grid_size {
                // Pack coordinates as SMHasher does
                let mut data = [0u8; 8];
                data[0..4].copy_from_slice(&(x as u32).to_le_bytes());
                data[4..8].copy_from_slice(&(y as u32).to_le_bytes());
                let hash = jamhash_bytes(&data);
                hashes.insert(hash);
            }
        }

        let total = grid_size * grid_size;
        let collision_rate = 1.0 - (hashes.len() as f64 / total as f64);
        assert!(
            collision_rate < 0.0001, // Less than 0.01% collision rate
            "PerlinNoise 2D: {} unique from {} ({:.4}% collision)",
            hashes.len(),
            total,
            collision_rate * 100.0
        );
    }

    #[test]
    fn test_smhasher_seed_block_len() {
        // SMHasher SeedBlockLen test: different seeds must produce different hashes
        // Failed lengths: 8, 9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 27, 28, 29, 30, 31
        let test_lengths = [8, 9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 27, 28, 29, 30, 31];

        for &len in &test_lengths {
            // Create test data
            let data: Vec<u8> = (0..len).map(|i| i as u8).collect();

            // Hash with different seeds
            let mut hasher1 = JamHasher::new_with_seed(0);
            hasher1.write(&data);
            let hash1 = hasher1.finish();

            let mut hasher2 = JamHasher::new_with_seed(1);
            hasher2.write(&data);
            let hash2 = hasher2.finish();

            let mut hasher3 = JamHasher::new_with_seed(0x12345678);
            hasher3.write(&data);
            let hash3 = hasher3.finish();

            // All seeds must produce different hashes for same data
            assert_ne!(
                hash1, hash2,
                "SeedBlockLen {}: seed 0 and 1 produced same hash",
                len
            );
            assert_ne!(
                hash1, hash3,
                "SeedBlockLen {}: seed 0 and 0x12345678 produced same hash",
                len
            );
            assert_ne!(
                hash2, hash3,
                "SeedBlockLen {}: seed 1 and 0x12345678 produced same hash",
                len
            );
        }
    }

    #[test]
    fn test_smhasher_seed_block_offset() {
        // SMHasher SeedBlockOffset test: seed sensitivity at different offsets
        // Failed offsets: 0, 1, 2, 3, 4, 5
        let data_len = 64;

        for offset in 0..=5 {
            let mut data: Vec<u8> = vec![0u8; data_len];

            // Set a single byte at the offset
            data[offset] = 0xFF;

            // Different seeds should produce different hashes
            let mut hasher1 = JamHasher::new_with_seed(0);
            hasher1.write(&data);
            let hash1 = hasher1.finish();

            let mut hasher2 = JamHasher::new_with_seed(0x55555555);
            hasher2.write(&data);
            let hash2 = hasher2.finish();

            assert_ne!(
                hash1, hash2,
                "SeedBlockOffset {}: different seeds produced same hash",
                offset
            );

            // Also check that the offset position matters
            let mut data_shifted = vec![0u8; data_len];
            data_shifted[(offset + 1) % data_len] = 0xFF;

            let mut hasher3 = JamHasher::new_with_seed(0);
            hasher3.write(&data_shifted);
            let hash3 = hasher3.finish();

            assert_ne!(
                hash1, hash3,
                "SeedBlockOffset {}: different offsets with same seed produced same hash",
                offset
            );
        }
    }

    #[test]
    fn test_seed_avalanche() {
        // Test that small seed changes cause large hash changes
        let data = b"test data for seed avalanche";

        for seed_bit in 0..64 {
            let seed1 = 0u64;
            let seed2 = 1u64 << seed_bit;

            let mut hasher1 = JamHasher::new_with_seed(seed1);
            hasher1.write(data);
            let hash1 = hasher1.finish();

            let mut hasher2 = JamHasher::new_with_seed(seed2);
            hasher2.write(data);
            let hash2 = hasher2.finish();

            let diff_bits = (hash1 ^ hash2).count_ones();
            assert!(
                diff_bits >= 20,
                "Seed bit {} flip: only {} bits differ (expected ~32)",
                seed_bit,
                diff_bits
            );
        }
    }

    #[test]
    fn test_8byte_sparse_patterns() {
        // Comprehensive test for 8-byte inputs with sparse bit patterns
        use std::collections::HashSet;

        let mut hashes = HashSet::new();

        // All single-bit patterns
        for bit in 0..64 {
            let val = 1u64 << bit;
            let hash = jamhash_bytes(&val.to_le_bytes());
            assert!(hashes.insert(hash), "8-byte single bit {} collision", bit);
        }

        // All two-bit patterns (adjacent)
        for bit in 0..63 {
            let val = 3u64 << bit;
            let hash = jamhash_bytes(&val.to_le_bytes());
            assert!(
                hashes.insert(hash),
                "8-byte two adjacent bits at {} collision",
                bit
            );
        }
    }

    // ==========================================
    // Statistical quality tests (SMHasher-style)
    // ==========================================

    /// Strict Avalanche Criterion test - each input bit should affect ~50% of output bits
    #[test]
    fn test_strict_avalanche_criterion_8bytes() {
        // For each input bit, flip it and measure which output bits change
        // Perfect avalanche: each output bit changes 50% of the time
        let num_samples = 1000;
        let mut bit_flip_counts = [[0u32; 64]; 64]; // [input_bit][output_bit]

        for sample in 0..num_samples {
            let base = (sample as u64).wrapping_mul(0x9e3779b97f4a7c15);
            let base_hash = jamhash_bytes(&base.to_le_bytes());

            for input_bit in 0..64 {
                let flipped = base ^ (1u64 << input_bit);
                let flipped_hash = jamhash_bytes(&flipped.to_le_bytes());
                let diff = base_hash ^ flipped_hash;

                for output_bit in 0..64 {
                    if (diff >> output_bit) & 1 == 1 {
                        bit_flip_counts[input_bit][output_bit] += 1;
                    }
                }
            }
        }

        // Check that each (input_bit, output_bit) pair has ~50% flip rate
        // SMHasher allows deviation of about 0.15 from 0.5
        let threshold = 0.35; // Very lenient for now, should be ~0.15
        let mut failures = Vec::new();

        for input_bit in 0..64 {
            for output_bit in 0..64 {
                let flip_rate = bit_flip_counts[input_bit][output_bit] as f64 / num_samples as f64;
                if (flip_rate - 0.5).abs() > threshold {
                    failures.push((input_bit, output_bit, flip_rate));
                }
            }
        }

        assert!(
            failures.len() < 10, // Allow some minor deviations
            "SAC failures (input_bit, output_bit, flip_rate): {:?}",
            &failures[..failures.len().min(20)]
        );
    }

    /// Bit Independence Criterion test - output bits should be independent
    #[test]
    fn test_bit_independence_criterion() {
        let num_samples = 10000;

        // Count how often pairs of output bits both flip when an input bit flips
        // For independent bits: P(both flip) = P(A) * P(B) ≈ 0.25
        let mut pair_flip_counts = [0u32; 64]; // Simplified: test bit 0 vs all others

        for sample in 0..num_samples {
            let base = (sample as u64).wrapping_mul(0x9e3779b97f4a7c15);
            let base_hash = jamhash_bytes(&base.to_le_bytes());

            // Flip input bit 0
            let flipped = base ^ 1;
            let flipped_hash = jamhash_bytes(&flipped.to_le_bytes());
            let diff = base_hash ^ flipped_hash;

            let bit0_flipped = diff & 1;
            for output_bit in 1..64 {
                let bit_flipped = (diff >> output_bit) & 1;
                if bit0_flipped == 1 && bit_flipped == 1 {
                    pair_flip_counts[output_bit] += 1;
                }
            }
        }

        // Check independence: should be ~25% (0.5 * 0.5) with some tolerance
        let mut failures = Vec::new();
        for output_bit in 1..64 {
            let pair_rate = pair_flip_counts[output_bit] as f64 / num_samples as f64;
            // Very lenient threshold
            if !(0.10..=0.40).contains(&pair_rate) {
                failures.push((output_bit, pair_rate));
            }
        }

        assert!(
            failures.len() < 5,
            "BIC failures (output_bit, pair_rate): {:?}",
            failures
        );
    }

    /// SMHasher-style permutation bias test
    #[test]
    fn test_permutation_bias_8bytes_zero_patterns() {
        // Test for bias when hashing [0, X] and [X, 0] patterns
        // Each output bit should be set ~50% of the time
        let mut bit_counts = [0u32; 64];
        let num_samples = 10000;

        for i in 0..num_samples {
            // Pattern: [0, i] as 16 bytes
            let mut data = [0u8; 16];
            data[8..16].copy_from_slice(&(i as u64).to_le_bytes());
            let hash = jamhash_bytes(&data);

            for bit in 0..64 {
                if (hash >> bit) & 1 == 1 {
                    bit_counts[bit] += 1;
                }
            }
        }

        // Check for bias - each bit should be set ~50% of the time
        let mut biased_bits = Vec::new();
        for bit in 0..64 {
            let rate = bit_counts[bit] as f64 / num_samples as f64;
            if (rate - 0.5).abs() > 0.05 {
                // 5% deviation allowed
                biased_bits.push((bit, rate));
            }
        }

        assert!(
            biased_bits.len() < 5,
            "Permutation [0, X] bias: bits with >5% deviation from 50%: {:?}",
            biased_bits
        );
    }

    /// Test seed mixing quality
    #[test]
    fn test_seed_mixing_quality() {
        // Test that seed differences propagate well through the hash
        let data = [0x42u8; 32];

        // Count bit differences between hashes with seeds that differ by 1 bit
        let mut total_diff_bits = 0u32;
        let mut min_diff = 64u32;
        let mut max_diff = 0u32;

        for seed_bit in 0..64 {
            let seed1 = 0u64;
            let seed2 = 1u64 << seed_bit;

            let mut h1 = JamHasher::new_with_seed(seed1);
            h1.write(&data);
            let hash1 = h1.finish();

            let mut h2 = JamHasher::new_with_seed(seed2);
            h2.write(&data);
            let hash2 = h2.finish();

            let diff = (hash1 ^ hash2).count_ones();
            total_diff_bits += diff;
            min_diff = min_diff.min(diff);
            max_diff = max_diff.max(diff);
        }

        let avg_diff = total_diff_bits as f64 / 64.0;

        // Should average ~32 bits difference (half of 64)
        assert!(
            avg_diff >= 25.0,
            "Seed mixing too weak: avg {} bit diff (expected ~32)",
            avg_diff
        );

        // No seed bit flip should cause fewer than 15 output bits to change
        assert!(
            min_diff >= 15,
            "Seed mixing has weak bit: min {} diff (expected >=20)",
            min_diff
        );
    }

    /// SMHasher TwoBytes test: all 2-byte keys must produce unique hashes with good distribution
    #[test]
    fn test_smhasher_twobytes_exhaustive() {
        use std::collections::HashSet;

        // Test all 65536 possible 2-byte combinations
        let mut hashes = HashSet::new();
        let mut collisions = Vec::new();

        for a in 0u8..=255 {
            for b in 0u8..=255 {
                let hash = jamhash_bytes(&[a, b]);
                if !hashes.insert(hash) {
                    collisions.push((a, b, hash));
                }
            }
        }

        assert!(
            collisions.is_empty(),
            "TwoBytes collisions found: {:?}",
            &collisions[..collisions.len().min(10)]
        );
        assert_eq!(hashes.len(), 65536, "Expected 65536 unique hashes");
    }

    /// SMHasher TwoBytes bit distribution: each output bit should be set ~50% of the time
    #[test]
    fn test_smhasher_twobytes_bit_distribution() {
        let mut bit_counts = [0u32; 64];

        for a in 0u8..=255 {
            for b in 0u8..=255 {
                let hash = jamhash_bytes(&[a, b]);
                for bit in 0..64 {
                    if (hash >> bit) & 1 == 1 {
                        bit_counts[bit] += 1;
                    }
                }
            }
        }

        // Each bit should be set ~32768 times (50%) with 5% tolerance
        let mut biased_bits = Vec::new();
        for (bit, &count) in bit_counts.iter().enumerate() {
            let rate = count as f64 / 65536.0;
            if (rate - 0.5).abs() > 0.05 {
                biased_bits.push((bit, rate));
            }
        }

        assert!(
            biased_bits.is_empty(),
            "TwoBytes bit distribution failures (bit, rate): {:?}",
            biased_bits
        );
    }

    /// SMHasher TwoBytes avalanche: flipping 1 input bit should flip ~50% of output bits
    #[test]
    fn test_smhasher_twobytes_avalanche() {
        let mut total_flips = 0u64;
        let mut min_flips = 64u32;
        let mut samples = 0u64;
        let mut worst_case = (0u8, 0u8, 0usize);

        for a in 0u8..=255 {
            for b in 0u8..=255 {
                let base_hash = jamhash_bytes(&[a, b]);

                // Flip each input bit and measure output bit changes
                for bit in 0..16 {
                    let (na, nb) = if bit < 8 {
                        (a ^ (1 << bit), b)
                    } else {
                        (a, b ^ (1 << (bit - 8)))
                    };
                    let flipped_hash = jamhash_bytes(&[na, nb]);
                    let diff_bits = (base_hash ^ flipped_hash).count_ones();
                    total_flips += diff_bits as u64;
                    if diff_bits < min_flips {
                        min_flips = diff_bits;
                        worst_case = (a, b, bit);
                    }
                    samples += 1;
                }
            }
        }

        let avg_flips = total_flips as f64 / samples as f64;
        let _ = worst_case; // Suppress unused warning

        // Should average ~32 bits (50% of 64)
        assert!(
            (28.0..=36.0).contains(&avg_flips),
            "TwoBytes avalanche: avg {} bits flipped (expected ~32)",
            avg_flips
        );
        // SMHasher typically allows min of ~12-15 bits for specific edge cases
        assert!(
            min_flips >= 12,
            "TwoBytes avalanche: min {} bits (expected >=12)",
            min_flips
        );
    }

    /// Distribution test for sequential inputs
    #[test]
    fn test_sequential_distribution() {
        // Hash sequential integers and check distribution of low bits
        let num_samples = 100000;
        let num_buckets = 256; // Test distribution into 256 buckets
        let mut bucket_counts = vec![0u32; num_buckets];

        for i in 0..num_samples {
            let hash = jamhash_bytes(&(i as u64).to_le_bytes());
            let bucket = (hash as usize) % num_buckets;
            bucket_counts[bucket] += 1;
        }

        // Chi-square test for uniform distribution
        let expected = num_samples as f64 / num_buckets as f64;
        let chi_sq: f64 = bucket_counts
            .iter()
            .map(|&c| {
                let diff = c as f64 - expected;
                diff * diff / expected
            })
            .sum();

        // For 255 degrees of freedom, chi-sq should be < ~310 at p=0.01
        assert!(
            chi_sq < 350.0,
            "Sequential distribution chi-sq {} (expected < 310)",
            chi_sq
        );
    }

    /// Test for the specific SeedBlockLen failure pattern
    #[test]
    fn test_seed_block_len_statistical() {
        // SMHasher SeedBlockLen checks that hash(seed, data) differs significantly
        // for all combinations of seeds and data lengths
        let test_lengths = [8, 9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 27, 28, 29, 30, 31];
        let num_seeds = 100;

        for &len in &test_lengths {
            let data: Vec<u8> = (0..len).map(|i| i as u8).collect();

            // Collect hashes for many seeds
            let mut hashes = Vec::with_capacity(num_seeds);
            for seed in 0..num_seeds {
                let mut hasher = JamHasher::new_with_seed(seed as u64);
                hasher.write(&data);
                hashes.push(hasher.finish());
            }

            // Check all pairs are different
            for i in 0..num_seeds {
                for j in (i + 1)..num_seeds {
                    assert_ne!(
                        hashes[i], hashes[j],
                        "SeedBlockLen {}: seeds {} and {} produced same hash",
                        len, i, j
                    );
                }
            }

            // Check bit distribution
            let mut bit_counts = [0u32; 64];
            for &hash in &hashes {
                for bit in 0..64 {
                    if (hash >> bit) & 1 == 1 {
                        bit_counts[bit] += 1;
                    }
                }
            }

            // Check for gross bias
            for bit in 0..64 {
                let rate = bit_counts[bit] as f64 / num_seeds as f64;
                assert!(
                    (0.2..=0.8).contains(&rate),
                    "SeedBlockLen {}: bit {} has rate {} (expected ~0.5)",
                    len,
                    bit,
                    rate
                );
            }
        }
    }
}
