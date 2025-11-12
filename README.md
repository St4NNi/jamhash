# jamhash

A fast, non-cryptographic hash function designed for genomics and general-purpose hashing.

## Features

- **`jamhash_u64`**: Optimized for hashing genomic k-mers encoded as u64 (2-bit encoding)
- **`JamHasher`**: Streaming hasher implementing `std::hash::Hasher` trait
- **Dual-path accumulation**: Strong avalanche properties with minimal overhead
- **Predictable**: `hash(0) == 0` for consistent behavior
- **Tested**: Low collision rates with canonical k-mers up to k=21

## Genomics Use Cases

The `jamhash_u64` function is specifically designed for bioinformatics workflows:

- **MinHash sketching** of DNA/RNA sequences
- **K-mer counting** and frequency analysis
- **Sequence similarity** estimation and comparison

```rust
use jamhash::jamhash_u64;

// Hash a k-mer encoded as u64 (2-bit encoding: A=00, C=01, G=10, T=11)
let kmer: u64 = 0x00011011;
let hash = jamhash_u64(kmer);
```

## General Purpose Hashing

```rust
use jamhash::{jamhash_bytes, JamHasher};
use std::hash::Hasher;

// Hash arbitrary bytes
let hash = jamhash_bytes(b"hello world");

// Streaming hasher
let mut hasher = JamHasher::new();
hasher.write(b"hello ");
hasher.write(b"world");
let hash = hasher.finish();
```

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
jamhash = "0.1"
```

## Inspirations

jamhash draws inspiration from:
- [rapidhash](https://github.com/Nicoshev/rapidhash) - Fast hashing with folded multiply
- [foldhash](https://github.com/orlp/foldhash) - Efficient folding techniques
- [ahash](https://github.com/tkaitchuck/aHash) - High-quality non-cryptographic hashing

## License

MIT OR Apache-2.0
