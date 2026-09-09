use std::fs;
use std::path::PathBuf;
use std::process::Command;

use tempfile::tempdir;

#[test]
fn test_binary_output_with_static_file() {
    // Input file
    let mut input = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    input.push("tests/data/test_input.tsv");

    // Reference output
    let mut expected = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    expected.push("tests/data/test_output.tsv");

    // Temporary directory automatically cleaned up when dropped
    let tmp_dir = tempdir().expect("Failed to create temp dir");
    let output = tmp_dir.path().join("output.tsv");

    // Run binary
    let result = Command::new(env!("CARGO_BIN_EXE_quickseg"))
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .output()
        .expect("Failed to run binary");

    assert!(
        result.status.success(),
        "Program exited with failure:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );

    // Read generated and expected outputs
    let generated =
        fs::read_to_string(&output).expect("Failed to read generated output");

    let expected =
        fs::read_to_string(&expected).expect("Failed to read expected output");

    assert_eq!(generated, expected);
}

#[test]
fn test_median_overflow_error() {
    let tmp = tempdir().unwrap();

    let input = tmp.path().join("input.bed");
    let output = tmp.path().join("output.tsv");

    fs::write(
        &input,
        "\
chr1\t0\t100\t500
chr1\t100\t200\t500
chr1\t200\t300\t500
chr1\t300\t400\t500
chr1\t400\t500\t500
",
    )
    .unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_quickseg"))
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .arg("--median")
        .arg("100") // deliberately too small
        .output()
        .unwrap();

    assert!(
        !result.status.success(),
        "Expected program to fail due to median overflow"
    );

    let stderr = String::from_utf8_lossy(&result.stderr);

    assert!(
        stderr.contains("Increase the value of --median."),
        "Unexpected error message:\n{}",
        stderr
    );
}
