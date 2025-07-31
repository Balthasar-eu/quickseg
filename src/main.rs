use clap::Parser;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

/// Simple program to count occurrences of values in the 4th column
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Input filename
    filename: String,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    // Initialize the helper array of 1000 elements to 0
    let mut helper = [0u32; 1000];

    // Open the file
    let file = File::open(&args.filename)?;
    let reader = io::BufReader::new(file);

    // Read line by line
    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result?;

        let columns: Vec<&str> = line.split('\t').collect();

        if columns.len() < 4 {
            eprintln!("Warning: Line {} has less than 4 columns, skipping.", line_num + 1);
            continue;
        }

        // Try to parse the 4th column as a float and convert to integer index
        match columns[3].trim().parse::<f64>() {
            Ok(num) => {
                let index = num as usize;
                if index < helper.len() {
                    helper[index] += 1;
                } else {
                    eprintln!("Warning: Index {} out of bounds at line {}, skipping.", index, line_num + 1);
                }
            }
            Err(e) => {
                eprintln!("Error parsing float in line {}: {}", line_num + 1, e);
            }
        }
    }

    // Print the resulting helper array (optional)
    for (i, count) in helper.iter().enumerate() {
        if *count > 0 {
            println!("Index {}: {}", i, count);
        }
    }

    Ok(())
}
