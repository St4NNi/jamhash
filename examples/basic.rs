use jamhash::{JamHasher, jamhash_bytes, jamhash_u64};
use std::hash::Hasher;

fn main() {
    println!("=== jamhash Examples ===\n");

    // Example 1: Hash a single u64 (optimized path)
    println!("1. Hashing single u64 values:");
    let kmer = 0x1234567890abcdef;
    let hash = jamhash_u64(kmer);
    println!("   jamhash_u64(0x{:016x}) = 0x{:016x}", kmer, hash);
    println!(
        "   jamhash_u64(0) = 0x{:016x} (predictable zero)",
        jamhash_u64(0)
    );
    println!();

    // Example 2: Hash bytes (one-shot)
    println!("2. Hashing byte slices:");
    let data = b"hello world";
    let hash = jamhash_bytes(data);
    println!(
        "   jamhash_bytes({:?}) = 0x{:016x}",
        std::str::from_utf8(data).unwrap(),
        hash
    );
    println!();

    // Example 3: Streaming hasher
    println!("3. Streaming hasher (same result as one-shot):");
    let mut hasher = JamHasher::new();
    hasher.write(b"hello ");
    hasher.write(b"world");
    let hash_streaming = hasher.finish();
    println!("   Streaming hash = 0x{:016x}", hash_streaming);
    println!("   Matches one-shot: {}", hash == hash_streaming);
    println!();

    // Example 4: Different input lengths produce different hashes
    println!("4. Length sensitivity:");
    for len in 0..=8 {
        let data = vec![0u8; len];
        let hash = jamhash_bytes(&data);
        println!("   {} zero bytes: 0x{:016x}", len, hash);
    }
    println!();

    // Example 5: Use with HashMap (via Hash trait)
    println!("5. Using with std::collections::HashMap:");
    use std::collections::HashMap;
    use std::hash::BuildHasherDefault;

    let mut map: HashMap<String, i32, BuildHasherDefault<JamHasher>> =
        HashMap::with_hasher(BuildHasherDefault::default());

    map.insert("foo".to_string(), 1);
    map.insert("bar".to_string(), 2);
    map.insert("baz".to_string(), 3);

    println!("   Created HashMap with JamHasher");
    println!("   map['foo'] = {}", map.get("foo").unwrap());
    println!();

    // Example 6: Avalanche demonstration
    println!("6. Avalanche effect (bit difference):");
    let pairs = [
        (b"test" as &[u8], b"Test" as &[u8]),
        (b"hello", b"Hello"),
        (b"abc", b"abd"),
    ];

    for (a, b) in &pairs {
        let hash_a = jamhash_bytes(a);
        let hash_b = jamhash_bytes(b);
        let diff_bits = (hash_a ^ hash_b).count_ones();
        println!(
            "   {:?} vs {:?}: {} bits differ",
            std::str::from_utf8(a).unwrap(),
            std::str::from_utf8(b).unwrap(),
            diff_bits
        );
    }
}
