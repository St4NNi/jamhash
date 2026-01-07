[![Crates.io](https://img.shields.io/crates/v/jamhash.svg)](https://crates.io/crates/jamhash)
[![Documentation](https://docs.rs/jamhash/badge.svg)](https://docs.rs/jamhash)
[![License](https://img.shields.io/crates/l/jamhash.svg)](https://github.com/St4NNi/jamhash#license)

# jamhash

A fast, non-cryptographic hash function designed for genomics and general-purpose hashing.

## Features

- **`jamhash_u64`**: Optimized for hashing genomic k-mers encoded as u64 (2-bit encoding)
- **`jamhash_bytes`**: Size-optimized byte hashing with paths for short (≤16), medium (17-128), long (129-4095), and ultra-long (≥4KB) inputs
- **`JamHasher`**: Streaming hasher implementing `std::hash::Hasher`
- **Dual-path accumulation**: Strong avalanche properties with minimal overhead
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

// Streaming std::hash::Hasher implementation
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

## Quality Verification

jamhash passes all 188 tests in [SMHasher3](https://gitlab.com/fwojcik/smhasher3), the comprehensive hash function test suite:

```
Overall result: pass (188 / 188 passed)
```

Full test results are available in [smhasher3_jamhash.txt](smhasher3_jamhash.txt).

## Inspirations

jamhash draws inspiration from:
- [rapidhash](https://github.com/Nicoshev/rapidhash)
- [foldhash](https://github.com/orlp/foldhash)
- [ahash](https://github.com/tkaitchuck/aHash)

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.
