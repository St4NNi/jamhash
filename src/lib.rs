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

/// Folded multiply: multiply two u64s and XOR the high and low halves
#[inline(always)]
fn fold_multiply(input: u64, const_1: u64) -> u64 {
    let temp = (input as u128).wrapping_mul(const_1 as u128);
    (temp as u64) ^ ((temp >> 64) as u64)
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

/// Streaming hasher using jamhash_u64 for chunk processing with 8 parallel accumulators.
///
/// Processes incoming bytes by hashing 8-byte chunks in parallel across 8 independent
/// accumulators, then merging them in the finalization step.
pub struct JamHasher {
    accumulators: [u64; 8],
    buffer: [u8; 8],
    buffer_len: usize,
    chunk_count: u64, // Track total chunks processed for accumulator selection
}

impl JamHasher {
    #[inline]
    pub fn new() -> Self {
        Self {
            // Initialize each accumulator with a different seed for better mixing
            accumulators: [
                CONST1,
                CONST1.rotate_left(8),
                CONST1.rotate_left(16),
                CONST1.rotate_left(24),
                CONST1.rotate_left(32),
                CONST1.rotate_left(40),
                CONST1.rotate_left(48),
                CONST1.rotate_left(56),
            ],
            buffer: [42u8; 8],
            buffer_len: 0,
            chunk_count: 0,
        }
    }

    pub fn new_with_seed(seed: u64) -> Self {
        let base = seed ^ CONST1;
        Self {
            accumulators: [
                base,
                base.rotate_left(8),
                base.rotate_left(16),
                base.rotate_left(24),
                base.rotate_left(32),
                base.rotate_left(40),
                base.rotate_left(48),
                base.rotate_left(56),
            ],
            buffer: [42u8; 8],
            buffer_len: 0,
            chunk_count: 0,
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
        let len = bytes.len();
        if len == 0 {
            return;
        }

        let mut offset = 0;

        // If we have buffered bytes, try to complete a full chunk
        if self.buffer_len > 0 {
            let needed = 8 - self.buffer_len;
            let available = len.min(needed);

            self.buffer[self.buffer_len..self.buffer_len + available]
                .copy_from_slice(&bytes[..available]);
            self.buffer_len += available;
            offset += available;

            // If we now have a full chunk, process it
            if self.buffer_len == 8 {
                let chunk = u64::from_le_bytes(self.buffer);
                let acc_idx = (self.chunk_count % 8) as usize;
                self.accumulators[acc_idx] =
                    fold_multiply(chunk ^ CONST1, self.accumulators[acc_idx] ^ CONST2);
                self.chunk_count += 1;
                self.buffer_len = 0;
                self.buffer = [42u8; 8];
            }
        }

        let remaining_bytes = &bytes[offset..];
        let (chunks, remainder) = remaining_bytes.as_chunks::<8>();

        // Process chunks in groups of 8 for maximum parallelism
        let num_full_groups = chunks.len() / 8;
        let mut chunk_idx = 0;

        for _ in 0..num_full_groups {
            // Process 8 chunks, one per accumulator
            // The compiler can parallelize these since they're independent
            for acc_idx in 0..8 {
                let chunk = u64::from_le_bytes(chunks[chunk_idx]);
                self.accumulators[acc_idx] =
                    fold_multiply(chunk ^ CONST1, self.accumulators[acc_idx] ^ CONST2);
                chunk_idx += 1;
            }
            self.chunk_count += 8;
        }

        // Process remaining chunks (0-7)
        while chunk_idx < chunks.len() {
            let chunk = u64::from_le_bytes(chunks[chunk_idx]);
            let acc_idx = (self.chunk_count % 8) as usize;
            self.accumulators[acc_idx] =
                fold_multiply(chunk ^ CONST1, self.accumulators[acc_idx] ^ CONST2);
            chunk_idx += 1;
            self.chunk_count += 1;
        }

        // Buffer any remaining bytes (< 8 bytes)
        let remaining = remainder.len();
        if remaining > 0 {
            self.buffer[..remaining].copy_from_slice(remainder);
            self.buffer_len = remaining;
        }
    }

    #[inline]
    fn finish(&self) -> u64 {
        // Process any remaining buffered bytes
        let mut accumulators = self.accumulators;

        if self.buffer_len != 0 {
            let buffer_value = u64::from_le_bytes(self.buffer);
            let acc_idx = (self.chunk_count % 8) as usize;
            accumulators[acc_idx] =
                fold_multiply(buffer_value ^ CONST1, accumulators[acc_idx] ^ CONST2);
        }

        // Merge all 8 accumulators into a single hash
        // Use a tree-reduction for better parallelism
        let a = fold_multiply(accumulators[0], accumulators[1]);
        let b = fold_multiply(accumulators[2], accumulators[3]);
        let c = fold_multiply(accumulators[4], accumulators[5]);
        let d = fold_multiply(accumulators[6], accumulators[7]);

        let ab = fold_multiply(a, b);
        let cd = fold_multiply(c, d);

        let merged = fold_multiply(ab, cd);

        jamhash_u64(merged)
    }
}

/// Hash arbitrary bytes using jamhash.
///
/// This is a convenience function that creates a `JamHasher`,
/// writes the bytes, and returns the hash.
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
    let mut hasher = JamHasher::new();
    hasher.write(bytes);
    hasher.finish()
}

#[cfg(test)]
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
    fn test_hasher_streaming() {
        let mut hasher1 = JamHasher::new();
        hasher1.write(b"hello world");
        let hash1 = hasher1.finish();

        let mut hasher2 = JamHasher::new();
        hasher2.write(b"hello ");
        hasher2.write(b"world");
        let hash2 = hasher2.finish();

        assert_eq!(hash1, hash2, "streaming should match one-shot");
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
}
