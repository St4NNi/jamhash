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
    accumulator: u64,
}

impl JamHasher {
    #[inline]
    pub fn new() -> Self {
        Self { accumulator: 0 }
    }

    pub fn new_with_seed(seed: u64) -> Self {
        Self { accumulator: seed }
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
        if len <= 16 {
            self.accumulator = hash_bytes_short(bytes, self.accumulator);
        } else {
            unsafe {
                // SAFETY: we checked that the length is > 16 bytes.
                self.accumulator = hash_bytes_long(bytes, self.accumulator);
            }
        }
    }

    #[inline]
    fn finish(&self) -> u64 {
        jamhash_u64(self.accumulator)
    }
}

/// Hashes strings <= 16 bytes, has unspecified behavior when bytes.len() > 16.
#[inline(always)]
fn hash_bytes_short(bytes: &[u8], accumulator: u64) -> u64 {
    let len = bytes.len();
    let mut s0 = accumulator;
    let mut s1 = CONST1;
    // XOR the input into s0, s1, then multiply and fold.
    if len >= 8 {
        s0 ^= u64::from_ne_bytes(bytes[0..8].try_into().unwrap());
        s1 ^= u64::from_ne_bytes(bytes[len - 8..].try_into().unwrap());
    } else if len >= 4 {
        s0 ^= u32::from_ne_bytes(bytes[0..4].try_into().unwrap()) as u64;
        s1 ^= u32::from_ne_bytes(bytes[len - 4..].try_into().unwrap()) as u64;
    } else if len > 0 {
        let lo = bytes[0];
        let mid = bytes[len / 2];
        let hi = bytes[len - 1];
        s0 ^= lo as u64;
        s1 ^= ((hi as u64) << 8) | mid as u64;
    }
    fold_multiply(s0, s1)
}

/// Load 8 bytes into a u64 word at the given offset.
///
/// # Safety
/// You must ensure that offset + 8 <= bytes.len().
#[inline(always)]
unsafe fn load(bytes: &[u8], offset: usize) -> u64 {
    // In most (but not all) cases this unsafe code is not necessary to avoid
    // the bounds checks in the below code, but the register allocation became
    // worse if I replaced those calls which could be replaced with safe code.
    unsafe { bytes.as_ptr().add(offset).cast::<u64>().read_unaligned() }
}

/// Hashes strings > 16 bytes.
///
/// # Safety
/// v.len() must be > 16 bytes.
#[cold]
#[inline(never)]
unsafe fn hash_bytes_long(mut v: &[u8], accumulator: u64) -> u64 {
    let mut s0 = accumulator;
    let mut s1 = s0.wrapping_add(CONST1);

    if v.len() > 128 {
        let mut s2 = s0.wrapping_add(CONST2);
        let mut s3 = s0.wrapping_add(CONST1);

        if v.len() > 256 {
            let mut s4 = s0.wrapping_add(CONST2);
            let mut s5 = s0.wrapping_add(CONST1);
            loop {
                unsafe {
                    // SAFETY: we checked the length is > 256, we index at most v[..96].
                    s0 = fold_multiply(load(v, 0) ^ s0, load(v, 48) ^ CONST1);
                    s1 = fold_multiply(load(v, 8) ^ s1, load(v, 56) ^ CONST1);
                    s2 = fold_multiply(load(v, 16) ^ s2, load(v, 64) ^ CONST2);
                    s3 = fold_multiply(load(v, 24) ^ s3, load(v, 72) ^ CONST1);
                    s4 = fold_multiply(load(v, 32) ^ s4, load(v, 80) ^ CONST2);
                    s5 = fold_multiply(load(v, 40) ^ s5, load(v, 88) ^ CONST1);
                }
                v = &v[96..];
                if v.len() <= 256 {
                    break;
                }
            }
            s0 ^= s4;
            s1 ^= s5;
        }

        loop {
            unsafe {
                // SAFETY: we checked the length is > 128, we index at most v[..64].
                s0 = fold_multiply(load(v, 0) ^ s0, load(v, 32) ^ CONST1);
                s1 = fold_multiply(load(v, 8) ^ s1, load(v, 40) ^ CONST1);
                s2 = fold_multiply(load(v, 16) ^ s2, load(v, 48) ^ CONST1);
                s3 = fold_multiply(load(v, 24) ^ s3, load(v, 56) ^ CONST1);
            }
            v = &v[64..];
            if v.len() <= 128 {
                break;
            }
        }
        s0 ^= s2;
        s1 ^= s3;
    }

    let len = v.len();
    unsafe {
        // SAFETY: our precondition ensures our length is at least 16, and the
        // above loops do not reduce the length under that. This protects our
        // first iteration of this loop, the further iterations are protected
        // directly by the checks on len.
        s0 = fold_multiply(load(v, 0) ^ s0, load(v, len - 16) ^ CONST1);
        s1 = fold_multiply(load(v, 8) ^ s1, load(v, len - 8) ^ CONST1);
        if len >= 32 {
            s0 = fold_multiply(load(v, 16) ^ s0, load(v, len - 32) ^ CONST1);
            s1 = fold_multiply(load(v, 24) ^ s1, load(v, len - 24) ^ CONST1);
            if len >= 64 {
                s0 = fold_multiply(load(v, 32) ^ s0, load(v, len - 48) ^ CONST2);
                s1 = fold_multiply(load(v, 40) ^ s1, load(v, len - 40) ^ CONST2);
                if len >= 96 {
                    s0 = fold_multiply(load(v, 48) ^ s0, load(v, len - 64) ^ CONST1);
                    s1 = fold_multiply(load(v, 56) ^ s1, load(v, len - 56) ^ CONST1);
                }
            }
        }
    }
    s0 ^ s1
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
