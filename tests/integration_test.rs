use std::path::PathBuf;
use std::process::Command;

#[test]
fn test_binary_output_with_static_file() {
    // Path to test file
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/data/test_input.tsv");

    // Run the compiled binary using `cargo run --bin <name>`
    let output = Command::new(env!("CARGO_BIN_EXE_quickseg")) // replace with your actual binary name
        .arg(&path)
        .output()
        .expect("Failed to run binary");

    assert!(
        output.status.success(),
        "Program exited with failure: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8_lossy(&output.stdout);

    // Check for expected output lines
    assert!(stdout.contains("Index 0: 1"));
    assert!(stdout.contains("Index 3: 4"));
    assert!(stdout.contains("Index 5: 1"));
    assert!(stdout.contains("Index 7: 3"));
    assert!(stdout.contains("Index 999: 1"));
}
