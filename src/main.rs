use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Read, Write, BufWriter, BufReader};
use std::path::PathBuf;
use flate2::read::MultiGzDecoder;

macro_rules! debug_println {
    ($($arg:tt)*) => (if ::std::cfg!(debug_assertions) { ::std::println!($($arg)*); })
}

/// Segment read counts into copy number segments
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Input file in TSV format. Needs to be in bed format with the fourth column as counts. You can use the output of mosdepth for this.
    #[arg(short, long, value_name = "FILE", value_parser = clap::value_parser!(PathBuf))]
    input: PathBuf,

    /// Output file to write segments
    #[arg(short, long, value_name = "FILE", value_parser = clap::value_parser!(PathBuf))]
    output: PathBuf,

    /// Optional output file for normalized counts that are used for segmenting. Useful for plotting or troubleshooting.
    #[arg(long, value_name = "FILE", value_parser = clap::value_parser!(PathBuf))]
    normalout: Option<PathBuf>,

    /// Optional normal sample input file
    #[arg(long, value_name = "FILE", value_parser = clap::value_parser!(PathBuf))]
    normal: Option<PathBuf>,

    /// Expected median. This needs to be higher than the median value of the raw counts. Increase if you use panel or amplicon sequencing.
    #[arg(long, default_value_t = 1000)]
    median: usize,

    /// Penalty parameter for segmentation. Lower values -> more segments. Higher values -> less segments
    #[arg(long, default_value_t = 10.0)]
    penalty: f64,

    /// Optional exclusion file for masking. Not implemented yet :(
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
                None, // Some(&PathBuf::from("normal.tsv"))
                true,
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
        normal_result,
        args.normalout.as_ref(),
        false,
    ).unwrap();


    let _ = write_to_tsv(&output, args.output);

    Ok(())
}


fn write_to_tsv(rows: &[TableRow], path: PathBuf) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(writer, "chr\tstart\tend\tvalue")?;

    // Write each row
    for row in rows {
        writeln!(writer, "{}\t{}\t{}\t{:.2}", row.chr, row.start, row.end, row.value)?;
    }

    Ok(())
}



fn open_reader(path: &PathBuf) -> std::io::Result<BufReader<Box<dyn Read>>> {
    let file = File::open(path)?;

    let reader: Box<dyn Read> = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    Ok(BufReader::new(reader))
    // Ok(BufReader::with_capacity(1024 * 1024, reader))
}



fn segment_file(
    input: &PathBuf,
    est_median: usize,
    penalty: f64,
    normal_segments: Option<Vec<TableRow>>,
    valuefile: Option<&PathBuf>,
    normal: bool,
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

    let mut reader  = open_reader(&input)?;

    // Read first and second line and save their byte offsets
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;
    let mut second_line = String::new();
    reader.read_line(&mut second_line)?;

    // Determine if first line is a header by checking if col2 is numeric
    let has_header = {
        let col2 = first_line.split('\t').nth(1).unwrap_or("").trim();
        col2.parse::<i32>().is_err()
    };

    // Determine bin size and starting chromosome
    let header_line = if has_header { &second_line } else { &first_line };
    let cols: Vec<&str> = header_line.split('\t').collect();

    let mut prev_chr = cols.get(0).unwrap_or(&"").trim().to_string();
    let start = cols.get(1).unwrap_or(&"").trim().parse::<usize>().unwrap_or(0);
    let end = cols.get(2).unwrap_or(&"").trim().parse::<usize>().unwrap_or(0);
    let bin_size = if start == 0 { Some(end) } else { None };

    let mut chrom_data: Vec<(String, (Vec<u32>, Vec<u32>, Vec<u32>))> = Vec::new();

    // Create and preallocate vectors
    let mut starts = Vec::new();
    let mut ends   = Vec::new();
    let mut values = Vec::new();
    
    let mut median_helper = vec![0u32; est_median];
    let mut element_count = 0;

    if let (Some(bs), Some(&chr_len)) = (bin_size, chrom_sizes.get(prev_chr.as_str())) {
        let bins = (chr_len + bs - 1) / bs;
        starts.reserve(bins);
        ends.reserve(bins);
        values.reserve(bins);
    }

    let mut result_table: Vec<TableRow> = Vec::new();

    let mut reader  = open_reader(&input)?;
    if has_header {
        let mut dummy = String::new();
        reader.read_line(&mut dummy)?; // consumes the first line
    }

    let lines = reader.lines();

    // Process data lines
    for (line_num, line_result) in lines.enumerate() {
        let line = line_result?;

        // TODO: This could possibly be better, but different things I tried either made performance worse or did nothing.
        // read_line instead of lines ~5 times slower, for whatever reason.
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
            debug_println!("Chromosome changed: {} → {}", prev_chr, chr);

            // Create new vectors for the new chromosome
            if let (Some(bs), Some(&chr_len)) = (bin_size, chrom_sizes.get(chr)) {
                let bins = (chr_len + bs - 1) / bs;
                starts = Vec::with_capacity(bins);
                ends = Vec::with_capacity(bins);
                values = Vec::with_capacity(bins);
                debug_println!("Preallocated {} bins for {} (bin size {})", bins, chr, bs);
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
                if !["chrX", "X", "chrY", "Y"].contains(&chr) {
                    element_count += 1;
                    if index < est_median {
                        median_helper[index] += 1;
                    } else {
                        debug_println!("Value {} larger than estimated median {} at line {}", index, est_median, line_num + 1);
                    }
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
    debug_println!("Median {}", median);

    // --- Normalize and pass to x() ---
    let mut loc_idx = 0;


    let mut writerx = match valuefile {
        Some(valuefile) => Some(BufWriter::new(File::create(valuefile)?)),
        None => None,
    };


    for (chr, (starts, ends, values)) in chrom_data.iter_mut() {

        median_helper.fill(0);
        let mut overflow_values = Vec::new();
        let mut newi = 0;

        //normal_segments[loc_idx].value
        for i in 0..values.len() {

            let raw_val = values[i] * 100;
            let normalized = match normal_segments {
                Some(ref normal_segments) => {
                    if *chr != normal_segments[loc_idx].chr {
                        loc_idx += 1;
                    }
                    if starts[i] >= normal_segments[loc_idx].end {
                        loc_idx += 1;
                    }
                    let normal_val = normal_segments[loc_idx].value;
                    if normal_val > 0.33 && normal_val < 3.00 {
                        Some( (raw_val * 100) / median / ((normal_val * 100.0) as u32) )
                    } else {
                        None
                    }
                }
                None => Some(raw_val / median),
            };


            if let Some(val) = normalized {
                if (val as usize) < est_median {
                    median_helper[val as usize] = 1;
                } else {
                    overflow_values.push(val as u32);
                    debug_println!("Value overflow {}", val);
                }

                if i != newi {
                    starts[newi] = starts[i];
                    ends[newi] = ends[i];
                }
                values[newi] = val;
                newi += 1;
            }
        }

    starts.truncate(newi);
    ends.truncate(newi);
    values.truncate(newi);

    // Write each row
    match writerx {
        Some(ref mut writerx) => { for i in 0..newi {writeln!(writerx, "{}\t{}\t{}\t{}", chr, starts[i], ends[i], values[i])?;}
        writerx.flush()?;
        }
        None => debug_println!("No output"),
    }

    let norm_values: Vec<f64> = values.iter().map(|&v| v as f64 / 100.0).collect();

    // Extract unique values (already "sorted" by index)
    let mut seg_values: Vec<f64> = median_helper
        .iter()
        .enumerate()
        .filter_map(|(i, &flag)| if flag == 1 { Some(i as f64 / 100.0) } else { None })
        .collect();

    overflow_values.sort_unstable();
    overflow_values.dedup();
    seg_values.extend(overflow_values.iter().map(|&v| v as f64 / 100.0));

    debug_println!("Segmentation {}", chr);

    let n = norm_values.len();
    let s = seg_values.len();

    if n == 0 {
        eprintln!("Skipping chromosome {}: no values after filtering", chr);
        continue;
    }

    let mut out_index = vec![0; n];
    let mut out_values = vec![0.0f64; n];

    let mut score = vec![0.0f64; s];
    // TODO: this should probably be a bitarray
    let mut backbool = vec![false; s*n];

    let mut backidx = vec![0; n];
    let mut breakidx = vec![0; n];

    score[0] = (norm_values[0] - seg_values[0]).abs();
    let mut minscore = score[0];

    for j in 1..s {
        score[j] = (norm_values[0] - seg_values[j]).abs();

        if minscore > score[j] {
            minscore = score[j];
            backidx[0] = j;
        }
    }

    minscore += penalty;

    for i in 1..n {
        if minscore < score[0] {
            backbool[s*i] = true;
            score[0] = minscore;
        }
        score[0] +=  (norm_values[i] - seg_values[0]).abs();
        let mut jmin = score[0];

        for j in 1..s {
            if minscore < score[j] {
                backbool[s*i+j] = true;
                score[j] = minscore;
            }
            score[j] += (norm_values[i] - seg_values[j]).abs();

            if jmin > score[j] {
                jmin = score[j];
                backidx[i] = j;
            }
        }

        minscore = jmin + penalty;
    }

    let mut b = 1;
    breakidx[0] = n;
    let mut maxixtmp = backidx[n-1];
    for i in (1..n).rev()  {
        if backbool[i*s+maxixtmp] {
            maxixtmp = backidx[i-1];
            breakidx[b] = i;
            b += 1;
        }
    }

    out_index[0] = 0;
    out_values[0] = seg_values[backidx[breakidx[b-1]-1]];

    for i in (0..b-1).rev() {
        out_index[b-(i+1)] = breakidx[i+1];
        out_values[b-(i+1)] = seg_values[backidx[breakidx[i]-1]];
    }

    out_index.truncate(b as usize);
    out_values.truncate(b as usize);
    
    for i in 0..out_index.len() {
        let start = starts[out_index[i] as usize];
        let end = if i + 1 < out_index.len() {
            starts[out_index[i + 1] as usize] - 1
        } else {
            *ends.last().unwrap()
        };
        let value = if normal && (chr == "chrX" || chr == "X" || chr == "chrY" || chr == "Y") {
             2.0 * out_values[i]
        } else {
            out_values[i]
        };

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

