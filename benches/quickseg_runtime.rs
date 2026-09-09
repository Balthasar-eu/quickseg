use criterion::{Criterion, criterion_group, criterion_main};
use std::process::Command;
use tempfile::tempdir;
use rand::{SeedableRng, rngs::StdRng};
use rand_distr::{Binomial, Distribution};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

fn generate_bed(path: &Path, multi: usize) {
    let mut rng = StdRng::seed_from_u64(42);

    let blocks = [
        (200, 40 * multi),
        (300, 10 * multi),
        (200, 20 * multi),
        (100, 15 * multi),
        (200, 15 * multi),
    ];

    let mut writer = BufWriter::new(File::create(path).unwrap());

    let mut pos = 0u64;

    for (n, len) in blocks {
        let dist = Binomial::new(n, 0.5).unwrap();

        for _ in 0..len {
            let value = dist.sample(&mut rng);

            writeln!(
                writer,
                "chr1\t{}\t{}\t{}",
                pos,
                pos + 100,
                value
            )
            .unwrap();

            pos += 100;
        }
    }
}


fn benchmark_quickseg(c: &mut Criterion) {
    let tmp = tempdir().unwrap();

    let input = tmp.path().join("input.bed");

    generate_bed(&input, 10_000);

    c.bench_function("quickseg_large", |b| {
        b.iter(|| {
            let output = tmp.path().join("out.seg");

            let result = Command::new(env!("CARGO_BIN_EXE_quickseg"))
                .arg("--input")
                .arg(&input)
                .arg("--output")
                .arg(&output)
                .output()
                .unwrap();

            assert!(result.status.success());

            let _ = std::fs::remove_file(&output);
        })
    });
}

criterion_group!(benches, benchmark_quickseg);
criterion_main!(benches);
