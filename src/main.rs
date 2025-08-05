use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Write, BufWriter, BufReader, Seek, SeekFrom};
use std::path::PathBuf;

/// Program to process BED-like file and count values by position
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Input file in TSV format
    #[arg(short, long, value_name = "FILE", value_parser = clap::value_parser!(PathBuf))]
    input: PathBuf,

    /// Output file to write segments
    #[arg(short, long, value_name = "FILE", value_parser = clap::value_parser!(PathBuf))]
    output: PathBuf,

    /// Optional normal sample input file
    #[arg(long, value_name = "FILE", value_parser = clap::value_parser!(PathBuf))]
    normal: Option<PathBuf>,

    /// Expected median (default: 1000)
    #[arg(long, default_value_t = 1000)]
    median: u32,

    /// Penalty parameter (default: 10.0)
    #[arg(long, default_value_t = 10.0)]
    penalty: f64,

    /// Optional exclusion file for masking
    #[arg(long, value_name = "FILE", value_parser = clap::value_parser!(PathBuf))]
    exclude: Option<PathBuf>,
}

#[derive(Debug, Clone)]
struct TableRow {
    chr: String,
    start: u32,
    end: u32,
    value: f64,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    let normal_result: Option<Vec<TableRow>> = if let Some(normal_path) = &args.normal {
        if normal_path.is_file() {
            Some(segment_file(
                normal_path,
                args.median,
                args.penalty,
                None,
                None,
            )?)
        } else {
            eprintln!("Provided --normal path is not a valid file: {:?}", normal_path);
            None
        }
    } else {
        None
    };
    
    let output = segment_file(
        &args.input,
        args.median,
        args.penalty,
        None, // mask: optional, to be handled inside if needed,
        None
    ).unwrap();


    write_to_tsv(&output, args.output);


    
    Ok(())
}


fn write_to_tsv(rows: &[TableRow], path: PathBuf) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(writer, "chr\tstart\tend\tvalue")?;

    // Write each row
    for row in rows {
        writeln!(writer, "{}\t{}\t{}\t{}", row.chr, row.start, row.end, row.value)?;
    }

    Ok(())
}


fn segment_file(
    input: &PathBuf,
    est_median: u32,
    penalty: f64,
    location_factors: Option<&Vec<u32>>,
    mask: Option<&Vec<bool>>,
) -> io::Result<Vec<TableRow>> {
    
    let chrom_sizes: HashMap<&str, usize> = [
        ("chr1", 249_000_000), ("1", 249_000_000),
        ("chr2", 242_200_000), ("2", 242_200_000),
        ("chr3", 198_300_000), ("3", 198_300_000),
        ("chr4", 190_300_000), ("4", 190_300_000),
        ("chr5", 181_600_000), ("5", 181_600_000),
        ("chr6", 170_900_000), ("6", 170_900_000),
        ("chr7", 159_400_000), ("7", 159_400_000),
        ("chr8", 145_200_000), ("8", 145_200_000),
        ("chr9", 138_400_000), ("9", 138_400_000),
        ("chr10", 133_800_000), ("10", 133_800_000),
        ("chr11", 135_100_000), ("11", 135_100_000),
        ("chr12", 133_300_000), ("12", 133_300_000),
        ("chr13", 114_400_000), ("13", 114_400_000),
        ("chr14", 107_100_000), ("14", 107_100_000),
        ("chr15", 102_000_000), ("15", 102_000_000),
        ("chr16", 90_400_000), ("16", 90_400_000),
        ("chr17", 83_300_000), ("17", 83_300_000),
        ("chr18", 80_400_000), ("18", 80_400_000),
        ("chr19", 58_700_000), ("19", 58_700_000),
        ("chr20", 64_500_000), ("20", 64_500_000),
        ("chr21", 46_800_000), ("21", 46_800_000),
        ("chr22", 50_900_000), ("22", 50_900_000),
        ("chrX", 156_100_000), ("X", 156_100_000),
        ("chrY", 57_300_000), ("Y", 57_300_000),
    ].into_iter().collect();

    // Open file with seek support
    let mut file = File::open(input)?;
    let mut reader = BufReader::new(&file);

    // Read first and second line and save their byte offsets
    let first_pos = reader.stream_position()?;
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;
    let second_pos = reader.stream_position()?;
    let mut second_line = String::new();
    reader.read_line(&mut second_line)?;

    // Determine if first line is a header by checking if col2 is numeric
    let has_header = {
        let col2 = first_line.split('\t').nth(1).unwrap_or("").trim();
        col2.parse::<i32>().is_err()
    };

    // Reset file position
    let seek_to = if has_header { second_pos } else { first_pos };
    file.seek(SeekFrom::Start(seek_to))?;
    let reader = BufReader::new(file);
    let lines = reader.lines();

    // Determine bin size and starting chromosome
    let header_line = if has_header { &second_line } else { &first_line };
    let cols: Vec<&str> = header_line.split('\t').collect();

    let chr = cols.get(0).unwrap_or(&"").trim().to_string();
    let start = cols.get(1).unwrap_or(&"").trim().parse::<usize>().unwrap_or(0);
    let end = cols.get(2).unwrap_or(&"").trim().parse::<usize>().unwrap_or(0);
    let bin_size = if start == 0 { Some(end) } else { None };

    //let mut chrom_data: HashMap<String, (Vec<u32>, Vec<u32>, Vec<u32>)> = HashMap::new();
    let mut chrom_data: Vec<(String, (Vec<u32>, Vec<u32>, Vec<u32>))> = Vec::new();

    // Create and preallocate vectors
    let mut starts = Vec::new();
    let mut ends   = Vec::new();
    let mut values = Vec::new();
    
    let mut median_helper = [0u32; 1000];
    let mut element_count = 0;

    if let (Some(bs), Some(&chr_len)) = (bin_size, chrom_sizes.get(chr.as_str())) {
        let bins = (chr_len + bs - 1) / bs;
        starts.reserve(bins);
        ends.reserve(bins);
        values.reserve(bins);
    }

    let mut prev_chr = chr;
    
    let mut result_table: Vec<TableRow> = Vec::new();

    // Process data lines
    for (line_num, line_result) in lines.enumerate() {
        let line = line_result?;
        let columns: Vec<&str> = line.split('\t').collect();

        if columns.len() < 4 {
            eprintln!("Skipping line {}: not enough columns", line_num + 1);
            continue;
        }

        let chr = columns[0].trim();

        // Detect chromosome change
        if chr != prev_chr {
            // Save previous chromosome data

            if !chrom_sizes.contains_key(chr) {
            continue; // Skip this iteration if key is not found.
            }

            chrom_data.push((prev_chr.clone(), (starts, ends, values)));
            println!("Chromosome changed: {} → {}", prev_chr, chr);

            // Create new vectors for the new chromosome
            if let (Some(bs), Some(&chr_len)) = (bin_size, chrom_sizes.get(chr)) {
                let bins = (chr_len + bs - 1) / bs;
                starts = Vec::with_capacity(bins);
                ends = Vec::with_capacity(bins);
                values = Vec::with_capacity(bins);
                println!("Preallocated {} bins for {} (bin size {})", bins, chr, bs);
            } else {
                starts = Vec::new();
                ends = Vec::new();
                values = Vec::new();
            }

            prev_chr = chr.to_string();
        }

        let start = columns[1].trim().parse::<u32>();
        let end = columns[2].trim().parse::<u32>();
        let value = columns[3].trim().split('.').next().unwrap().parse::<u32>();

        match (start, end, value.clone()) {
            (Ok(s), Ok(e), Ok(v)) => {
                let index = v as usize;
                starts.push(s);
                ends.push(e);
                values.push(v);
                element_count += 1;
                if index < 1000 {
                    median_helper[index] += 1;
                } else {
                    eprintln!("Value index {} out of range at line {}", index, line_num + 1);
                }
            }
            _ => {
                eprintln!("Invalid number at line {}", line_num + 1);
            }
        }
    }

    chrom_data.push((prev_chr, (starts, ends, values)));

    // --- Median calculation ---
    let mut cumulative = 0;
    let mut median = 0;

    for (i, count) in median_helper.iter().enumerate() {
        cumulative += count;
        if cumulative as f64 >= element_count as f64 / 2.0 {
            median = i as u32;
            break;
        }
    }

    // --- Normalize and pass to x() ---
    for (chr, (starts, ends, values)) in &chrom_data {

    // Constants
    const MAX_BUCKETS: usize = 1000;

    // Prepare data structures
    let mut masked_starts = Vec::with_capacity(values.len());
    let mut masked_ends = Vec::with_capacity(values.len());
    let mut norm_values = Vec::with_capacity(values.len());
    let mut bucket_flags = [0u8; MAX_BUCKETS]; // Used as bitmap for unique values
    let mut overflow_values = Vec::new();

    let mut loc_idx = 0;

    for i in 0..values.len() {
        let use_entry = mask.map_or(true, |m| m[i]); // If mask is None, include all entries

        if use_entry {
            let raw_val = values[i] * 100;
            let normalized = match location_factors {
                Some(factors) => raw_val / median / factors[loc_idx],
                None => raw_val / median,
            };

            if (normalized as usize) < MAX_BUCKETS {
                bucket_flags[normalized as usize] = 1;
            } else {
                overflow_values.push(normalized as usize);
            }

            masked_starts.push(starts[i]);
            masked_ends.push(ends[i]);
            norm_values.push(normalized as f64 / 100.0);

            loc_idx += 1;
        }
    }

    // Extract unique values (already "sorted" by index)
    let mut seg_values: Vec<f64> = bucket_flags
        .iter()
        .enumerate()
        .filter_map(|(i, &flag)| if flag == 1 { Some(i as f64 / 100.0) } else { None })
        .collect();

    overflow_values.sort_unstable();
    overflow_values.dedup();
    seg_values.extend(overflow_values.iter().map(|&v| v as f64 / 100.0));
        
    let mut out_index = vec![0; norm_values.len()];
    let mut out_values = vec![0.0f64; norm_values.len()];
    let outsize = unsafe {
        segment(norm_values.as_ptr(), norm_values.len(),  seg_values.as_ptr(), seg_values.len(), penalty, out_index.as_mut_ptr(), out_values.as_mut_ptr())
    };
    out_index.truncate(outsize as usize);
    out_values.truncate(outsize as usize);
    
    for i in 0..out_index.len() {
        let start = masked_starts[out_index[i] as usize];
        let end = if i + 1 < out_index.len() {
            masked_ends[out_index[i + 1] as usize]
        } else {
            *masked_ends.last().unwrap()
        };
        let value = out_values[i];

        result_table.push(TableRow {
            chr: chr.clone(),
            start,
            end,
            value,
        });
    }
        
    }
    Ok(result_table)
}


unsafe extern "C" {
    fn segment(
        val_values: *const f64,
        n: usize,
        seg_values: *const f64,
        s: usize,
        penalty: f64,
        out_index: *mut i32,
        out_values: *mut f64,
    ) -> i32;
}
