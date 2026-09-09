use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Read, Write, BufWriter, BufReader};
use std::path::PathBuf;
use flate2::read::MultiGzDecoder;
// use std::time::Instant;
use atoi::FromRadix10;

mod segmentation;
use segmentation::segment_forward;

mod bufreader;
use bufreader::CompactBufReader;

macro_rules! debug_println {
    ($($arg:tt)*) => (if ::std::cfg!(debug_assertions) { ::std::println!($($arg)*); })
}

#[inline]
fn get_bit(bits: &[u64], idx: usize) -> bool {
    (bits[idx / 64] & (1u64 << (idx % 64))) != 0
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
    )?;


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

fn open_reader_c(path: &PathBuf) -> std::io::Result<CompactBufReader<Box<dyn Read>>> {
    let file = File::open(path)?;

    let reader: Box<dyn Read> = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    //Ok(CompactBufReader::new(reader))
    Ok(CompactBufReader::with_capacity(reader, 8 * 64 * 1024))
}

fn segment_file(
    input: &PathBuf,
    est_median: usize,
    penalty: f64,
    normal_segments: Option<Vec<TableRow>>,
    valuefile: Option<&PathBuf>,
    normal: bool,
) -> io::Result<Vec<TableRow>> {
    // let now = Instant::now();

    // TODO: the hg19 chrom sizes are wrong! copy from here: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes
    let chrom_sizes: HashMap<&'static [u8], usize> = HashMap::from([
        (b"chr1" as &'static [u8], 249_000_000), (b"1" as &'static [u8], 249_000_000),
        (b"chr2" as &'static [u8], 242_200_000), (b"2" as &'static [u8], 242_200_000),
        (b"chr3" as &'static [u8], 198_300_000), (b"3" as &'static [u8], 198_300_000),
        (b"chr4" as &'static [u8], 190_300_000), (b"4" as &'static [u8], 190_300_000),
        (b"chr5" as &'static [u8], 181_600_000), (b"5" as &'static [u8], 181_600_000),
        (b"chr6" as &'static [u8], 170_900_000), (b"6" as &'static [u8], 170_900_000),
        (b"chr7" as &'static [u8], 159_400_000), (b"7" as &'static [u8], 159_400_000),
        (b"chr8" as &'static [u8], 145_200_000), (b"8" as &'static [u8], 145_200_000),
        (b"chr9" as &'static [u8], 138_400_000), (b"9" as &'static [u8], 138_400_000),
        (b"chr10" as &'static [u8], 133_800_000), (b"10" as &'static [u8], 133_800_000),
        (b"chr11" as &'static [u8], 135_100_000), (b"11" as &'static [u8], 135_100_000),
        (b"chr12" as &'static [u8], 133_300_000), (b"12" as &'static [u8], 133_300_000),
        (b"chr13" as &'static [u8], 114_400_000), (b"13" as &'static [u8], 114_400_000),
        (b"chr14" as &'static [u8], 107_100_000), (b"14" as &'static [u8], 107_100_000),
        (b"chr15" as &'static [u8], 102_000_000), (b"15" as &'static [u8], 102_000_000),
        (b"chr16" as &'static [u8], 90_400_000), (b"16" as &'static [u8], 90_400_000),
        (b"chr17" as &'static [u8], 83_300_000), (b"17" as &'static [u8], 83_300_000),
        (b"chr18" as &'static [u8], 80_400_000), (b"18" as &'static [u8], 80_400_000),
        (b"chr19" as &'static [u8], 58_700_000), (b"19" as &'static [u8], 58_700_000),
        (b"chr20" as &'static [u8], 64_500_000), (b"20" as &'static [u8], 64_500_000),
        (b"chr21" as &'static [u8], 46_800_000), (b"21" as &'static [u8], 46_800_000),
        (b"chr22" as &'static [u8], 50_900_000), (b"22" as &'static [u8], 50_900_000),
        (b"chrX" as &'static [u8], 156_100_000), (b"X" as &'static [u8], 156_100_000),
        (b"chrY" as &'static [u8], 57_300_000), (b"Y" as &'static [u8], 57_300_000),
    ]);

    let mut reader  = open_reader(&input)?;

    // Read first and second line and save their byte offsets
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;
    let mut second_line = String::new();
    reader.read_line(&mut second_line)?;

    // Determine if first line is a header by checking if col2 is numeric
    let mut skipline = {
        let col2 = first_line.split('\t').nth(1).unwrap_or("").trim();
        col2.parse::<i32>().is_err()
    };

    // Determine bin size and starting chromosome
    let header_line = if skipline { &second_line } else { &first_line };
    let cols: Vec<&str> = header_line.split('\t').collect();

    let binding = cols.get(0).unwrap_or(&"").trim().to_string(); // TODO: this looks pretty stupid, but it works.
    let mut prev_chr = binding.as_bytes().to_vec();  // TODO: this looks pretty stupid, but it works.
    let start = cols.get(1).unwrap_or(&"").trim().parse::<usize>().unwrap_or(0);
    let end = cols.get(2).unwrap_or(&"").trim().parse::<usize>().unwrap_or(0);
    let bin_size = if start == 0 { Some(end) } else { None };

    let mut chrom_data: Vec<(String, (Vec<u32>, Vec<u32>, Vec<u32>))> = Vec::with_capacity(24);

    // Create and preallocate vectors
    let mut starts = Vec::new();
    let mut ends   = Vec::new();
    let mut values = Vec::new();
    
    let mut median_helper = vec![0u32; est_median];
    let mut element_count = 0;

    if let (Some(bs), Some(&chr_len)) = (bin_size, chrom_sizes.get(prev_chr.as_slice())) {
        let bins = (chr_len + bs - 1) / bs;
        starts.reserve(bins);
        ends.reserve(bins);
        values.reserve(bins);
    }

    let mut result_table: Vec<TableRow> = Vec::new();

    let mut rdr = open_reader_c(&input)?;

    // println!("Elapsed start: {:.2?}", now.elapsed());

    let mut linecount = 0usize;
    let mut fieldcount = 0usize;
    let mut is_sex_chr = false;

    loop {
        if !rdr.ensure_delimiter()? {
            break;
        }

        while let Some(pos) = rdr.find_delimiter() {
            let cur = rdr.available();

            let field = &cur[..pos];
            let delim = cur[pos];

            if !skipline {
                match fieldcount {
                    0 => {
                        if field != prev_chr.as_slice() {
                            if !chrom_sizes.contains_key(field) {
                                skipline = true;
                            } else {
                                chrom_data.push((String::from_utf8(prev_chr).unwrap(), (starts, ends, values)));
                                // Create new vectors for the new chromosome
                                let mut bins = 1000;
                                if let (Some(bs), Some(&chr_len)) = (bin_size, chrom_sizes.get(field)) {
                                    bins = (chr_len + bs - 1) / bs;
                                }
                                starts = Vec::with_capacity(bins);
                                ends = Vec::with_capacity(bins);
                                values = Vec::with_capacity(bins);
                                prev_chr = field.to_vec();

                                is_sex_chr = matches!(field, b"chrX" | b"X" | b"chrY" | b"Y");
                            }
                        }

                    }
                    1 => {
                        let (s, _used) = u32::from_radix_10(field);
                        starts.push(s);
                    }
                    2 => {
                        let (e, _used) = u32::from_radix_10(field);
                        ends.push(e);
                    }
                    3 => {
                        let (v, _used) = u32::from_radix_10(field);
                        // println!("{},{}",v,_used);
                        let index = v as usize;
                        values.push(v);
                        if !is_sex_chr {
                            element_count += 1;
                            if index < est_median {
                                median_helper[index] += 1;
                            }
                            else {
                                debug_println!("Value {} larger than estimated median {} at line {}", index, est_median, linecount + 1);
                            }
                        }
                    }
                    _ => {println!("Error, {}. field!", fieldcount);}
                }
                fieldcount += 1;
            }

            rdr.consume_front(pos + 1);

            if delim == b'\n' {
                linecount += 1;
                fieldcount = 0;
                break;
            }
        }
    }

    // println!("Elapsed file read: {:.2?}", now.elapsed());

    chrom_data.push((String::from_utf8(prev_chr).unwrap(), (starts, ends, values)));

    // --- Median calculation ---
    let mut cumulative = 0;
    let mut median: Option<u32> = None;

    for (i, count) in median_helper.iter().enumerate() {
        cumulative += count;

        if cumulative as f64 >= element_count as f64 / 2.0 {
            median = Some(i as u32);
            break;
        }
    }

    let median = median.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Median exceeds configured limit (--median {}). \
                Increase the value of --median.",
                est_median
            ),
        )
    })?;

    debug_println!("Median {}", median);

    // --- Normalize and pass to x() ---
    let mut loc_idx = 0;


    let mut writerx = match valuefile {
        Some(valuefile) => Some(BufWriter::new(File::create(valuefile)?)),
        None => None,
    };

    // println!("Elapsed pre seg: {:.2?}", now.elapsed());


    for (chr, (starts, ends, values)) in chrom_data.iter_mut() {

        let normalized = normalize_chromosome(
            chr,
            starts,
            ends,
            values,
            normal_segments.as_ref(),
            &mut loc_idx,
            median,
            est_median,
            &mut writerx,
        )?;

        if normalized.norm_values.is_empty() {
            continue;
        }

        let (out_index, out_values) = segment_chromosome(
            &normalized.norm_values,
            &normalized.seg_values,
            penalty,
        );

        build_segments(
            &mut result_table,
            chr,
            starts,
            ends,
            &out_index,
            &out_values,
            normal,
        );
    }

    // println!("Elapsed post result: {:.2?}", now.elapsed());
    Ok(result_table)
}

struct NormalizedChromosome {
    seg_values: Vec<f64>,
    norm_values: Vec<f64>,
}

fn normalize_chromosome(
    chr: &str,
    starts: &mut Vec<u32>,
    ends: &mut Vec<u32>,
    values: &mut Vec<u32>,
    normal_segments: Option<&Vec<TableRow>>,
    loc_idx: &mut usize,
    median: u32,
    est_median: usize,
    writerx: &mut Option<BufWriter<File>>,
) -> io::Result<NormalizedChromosome> {

    let mut median_helper = vec![0u32; est_median];
    let mut overflow_values = Vec::new();
    let mut newi = 0;

    for i in 0..values.len() {
        let raw_val = values[i] * 100;

        let normalized = match normal_segments {
            Some(normal_segments) => {
                if chr != normal_segments[*loc_idx].chr {
                    *loc_idx += 1;
                }

                while starts[i] >= normal_segments[*loc_idx].end {
                    *loc_idx += 1;
                }

                let normal_val = normal_segments[*loc_idx].value;

                if normal_val > 0.33 && normal_val < 3.00 {
                    Some(
                        (raw_val * 100)
                            / median
                            / ((normal_val * 100.0) as u32),
                    )
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
                overflow_values.push(val);
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

    if let Some(writer) = writerx {
        for i in 0..newi {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                chr,
                starts[i],
                ends[i],
                values[i]
            )?;
        }
    }

    let norm_values =
        values.iter().map(|&v| v as f64 / 100.0).collect();

    let mut seg_values: Vec<f64> = median_helper
        .iter()
        .enumerate()
        .filter_map(|(i, &flag)| {
            if flag == 1 {
                Some(i as f64 / 100.0)
            } else {
                None
            }
        })
        .collect();

    overflow_values.sort_unstable();
    overflow_values.dedup();

    seg_values.extend(
        overflow_values
            .iter()
            .map(|&v| v as f64 / 100.0),
    );

    let padded_len = seg_values.len().next_multiple_of(8);
    seg_values.resize(padded_len, f64::INFINITY);

    Ok(NormalizedChromosome {
        seg_values,
        norm_values,
    })
}

fn segment_chromosome(
    norm_values: &[f64],
    seg_values: &[f64],
    penalty: f64,
) -> (Vec<usize>, Vec<f64>) {

    let n = norm_values.len();
    let s = seg_values.len();

    let mut out_index = vec![0; n];
    let mut out_values = vec![0.0; n];

    let num_bits = s * n;
    let mut backbool = vec![0u64; num_bits.div_ceil(64)];

    let mut backidx = vec![0usize; n];
    let mut breakidx = vec![0usize; n];

    segment_forward(
        norm_values,
        seg_values,
        penalty,
        &mut backbool,
        &mut backidx,
    );

    let mut b = 1;

    breakidx[0] = n;
    let mut maxixtmp = backidx[n - 1];

    for i in (1..n).rev() {
        if get_bit(&backbool, i * s + maxixtmp) {
            maxixtmp = backidx[i - 1];
            breakidx[b] = i;
            b += 1;
        }
    }

    out_index[0] = 0;
    out_values[0] = seg_values[backidx[breakidx[b - 1] - 1]];

    for i in (0..b - 1).rev() {
        out_index[b - (i + 1)] = breakidx[i + 1];
        out_values[b - (i + 1)] =
            seg_values[backidx[breakidx[i] - 1]];
    }

    out_index.truncate(b);
    out_values.truncate(b);

    (out_index, out_values)
}


fn build_segments(
    result_table: &mut Vec<TableRow>,
    chr: &str,
    starts: &[u32],
    ends: &[u32],
    out_index: &[usize],
    out_values: &[f64],
    normal: bool,
) {
    for i in 0..out_index.len() {
        let start = starts[out_index[i]];

        let end = if i + 1 < out_index.len() {
            starts[out_index[i + 1]] - 1
        } else {
            *ends.last().unwrap()
        };

        let value = if normal
            && matches!(chr, "chrX" | "X" | "chrY" | "Y")
        {
            2.0 * out_values[i]
        } else {
            out_values[i]
        };

        result_table.push(TableRow {
            chr: chr.to_string(),
            start,
            end,
            value,
        });
    }
}

