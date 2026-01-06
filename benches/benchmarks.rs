use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::{hash::BuildHasher, hint::black_box, time::Duration};

#[inline]
pub fn murmur3_old(value: u64) -> u64 {
    murmurhash3::murmurhash3_x64_128(&value.to_le_bytes(), 42).0
}

#[inline]
pub fn murmur3_new(value: u64) -> u64 {
    fastmurmur3::murmur3_x64_128(&value.to_le_bytes(), 42) as u64
}

#[inline]
pub fn xxhash3(value: u64) -> u64 {
    xxhash_rust::xxh3::xxh3_64(&value.to_le_bytes())
}

#[inline]
pub fn jamhash(value: u64) -> u64 {
    jamhash::jamhash_u64(value)
}

#[inline]
pub fn foldhash(value: u64) -> u64 {
    use std::hash::Hasher;
    let mut hasher = foldhash::quality::RandomState::default().build_hasher();
    hasher.write_u64(value);
    hasher.finish()
}

fn bench_hash_functions(c: &mut Criterion) {
    let mut group = c.benchmark_group("single_hash");
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    // Pre-generate enough values for the benchmark
    let mut values = (0..100000u64).cycle();

    group.bench_function("xxhash3", |b| {
        b.iter(|| xxhash3(black_box(values.next().unwrap())))
    });

    group.bench_function("murmur3_old", |b| {
        b.iter(|| murmur3_old(black_box(values.next().unwrap())))
    });

    group.bench_function("murmur3_new", |b| {
        b.iter(|| murmur3_new(black_box(values.next().unwrap())))
    });

    group.bench_function("jamhash", |b| {
        b.iter(|| jamhash(black_box(values.next().unwrap())))
    });

    group.bench_function("foldhash", |b| {
        b.iter(|| foldhash(black_box(values.next().unwrap())))
    });

    group.finish();
}

fn bench_bytes_hashing(c: &mut Criterion) {
    use std::hash::Hasher;

    // Test various sizes: short (≤16), medium (17-128), long (>128), ultra-long (>4KB)
    let sizes = [4, 8, 16, 32, 64, 128, 256, 1024, 4096, 8192, 16384];

    for &size in &sizes {
        let data: Vec<u8> = (0..size).map(|i| i as u8).collect();

        let mut group = c.benchmark_group(format!("bytes_{}", size));
        group.throughput(Throughput::Bytes(size as u64));
        group.warm_up_time(Duration::from_millis(300));
        group.measurement_time(Duration::from_secs(2));

        group.bench_function("jamhash", |b| {
            b.iter(|| jamhash::jamhash_bytes(black_box(&data)))
        });

        group.bench_function("xxhash3", |b| {
            b.iter(|| xxhash_rust::xxh3::xxh3_64(black_box(&data)))
        });

        group.bench_function("foldhash", |b| {
            b.iter(|| {
                let mut hasher = foldhash::quality::RandomState::default().build_hasher();
                hasher.write(black_box(&data));
                hasher.finish()
            })
        });

        group.finish();
    }
}

criterion_group!(benches, bench_hash_functions, bench_bytes_hashing);
criterion_main!(benches);
