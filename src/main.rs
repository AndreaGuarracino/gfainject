use anyhow::Result;
use roaring::RoaringBitmap;
use std::collections::BTreeMap;
use std::io::prelude::*;
use std::{io::BufReader, path::PathBuf};

use std::fs::File;
use gbam_tools::reader::parse_tmplt::ParsingTemplate;
use gbam_tools::reader::reader::Reader;
use gbam_tools::reader::records::Records;
use gbam_tools::Fields;

use clap::{Parser, ArgGroup};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(group(
    ArgGroup::new("input_mode")
        .required(true)
        .args(["bam", "paf", "gbam", "range"]),
))]
struct Args {
    /// Path to input GFA file
    #[arg(short, long, required = true)]
    gfa: PathBuf,

    /// Path to input BAM file
    #[arg(short, long)]
    bam: Option<PathBuf>,

    /// Path to input PAF file
    #[arg(short, long)]
    paf: Option<PathBuf>,

    /// Path to input GBAM file
    #[arg(short, long)]
    gbam: Option<PathBuf>,

    /// Range query in format "path_name:start-end"
    #[arg(short, long, value_parser = parse_range)]
    range: Option<(String, usize, usize)>,
}

/// Parse a range string in the format "path_name:start-end"
/// where path_name can contain ':' and '-' characters
fn parse_range(s: &str) -> Result<(String, usize, usize), String> {
    // Find the last colon in the string
    let last_colon_pos = s.rfind(':').ok_or("Range must include ':' separator".to_string())?;
    
    // Split into path_name and range parts
    let path_name = s[..last_colon_pos].to_string();
    let range_str = &s[last_colon_pos + 1..];

    // Split the range part by the last hyphen
    let range_parts: Vec<&str> = range_str.split('-').collect();
    if range_parts.len() != 2 {
        return Err("Range must be in format start-end".to_string());
    }

    // Parse start and end positions
    let start = range_parts[0].parse::<usize>()
        .map_err(|_| "Invalid start position".to_string())?;
    let end = range_parts[1].parse::<usize>()
        .map_err(|_| "Invalid end position".to_string())?;

    if end <= start {
        return Err("End position must be greater than start position".to_string());
    }

    Ok((path_name, start, end))
}

/// Represents a single step in a path through the graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct PathStep {
    /// Node ID in the path
    node: u32,
    /// Whether this node is traversed in reverse
    reverse: bool,
}

/// Index structure for efficient path queries
struct PathIndex {
    /// Range of segment IDs (min, max)
    segment_id_range: (usize, usize),
    /// Length of each segment
    segment_lens: Vec<usize>,
    /// Map from path names to path indices
    path_names: BTreeMap<String, usize>,
    /// Steps for each path
    path_steps: Vec<Vec<PathStep>>,
    /// Offset positions for each step in each path
    path_step_offsets: Vec<roaring::RoaringBitmap>,
}

struct PathStepRangeIter<'a> {
    path_id: usize,
    pos_range: std::ops::Range<u32>,
    // start_pos: usize,
    // end_pos: usize,
    steps: Box<dyn Iterator<Item = (usize, &'a PathStep)> + 'a>,
    // first_step_start_pos: u32,
    // last_step_end_pos: u32,
}

impl<'a> Iterator for PathStepRangeIter<'a> {
    type Item = (usize, &'a PathStep);

    fn next(&mut self) -> Option<Self::Item> {
        self.steps.next()
    }
}

impl PathIndex {
    fn path_step_range_iter<'a>(
        &'a self,
        path_name: &str,
        pos_range: std::ops::Range<u32>,
    ) -> Option<PathStepRangeIter<'a>> {
        let path_id = *self.path_names.get(path_name)?;
        let offsets = self.path_step_offsets.get(path_id)?;

        let start = pos_range.start;
        let end = pos_range.end;
        let start_rank = offsets.rank(start);
        let end_rank = offsets.rank(end);

        let steps = {
            let path_steps = self.path_steps.get(path_id)?;

            let skip = (start_rank as usize).checked_sub(1).unwrap_or_default();
            let take = end_rank as usize - skip;
            let iter = path_steps
                .iter()
                .skip(skip)
                .take(take)
                .enumerate()
                .map(move |(ix, step)| (skip + ix, step))
                .fuse();

            Box::new(iter) as Box<dyn Iterator<Item = _>>
        };

        Some(PathStepRangeIter {
            path_id,
            pos_range,
            steps,
            // first_step_start_pos,
            // last_step_end_pos,
        })
    }

    fn from_gfa(gfa_path: impl AsRef<std::path::Path>) -> Result<Self> {
        let gfa = std::fs::File::open(&gfa_path)?;
        let mut gfa_reader = BufReader::new(gfa);

        let mut line_buf = Vec::new();

        let mut seg_lens = Vec::new();

        let mut seg_id_range = (std::usize::MAX, 0usize);
        // dbg!();

        loop {
            line_buf.clear();

            let len = gfa_reader.read_until(0xA, &mut line_buf)?;
            if len == 0 {
                break;
            }

            let line = &line_buf[..len];
            let line_str = std::str::from_utf8(&line)?;
            // println!("{line_str}");

            if !matches!(line.first(), Some(b'S')) {
                continue;
            }

            let mut fields = line_str.split(|c| c == '\t');

            let Some((name, seq)) = fields.next().and_then(|_type| {
                let name = fields.next()?.trim();
                let seq = fields.next()?.trim();
                Some((name, seq))
            }) else {
                continue;
            };
            let seg_id = name.parse::<usize>()?;

            seg_id_range.0 = seg_id_range.0.min(seg_id);
            seg_id_range.1 = seg_id_range.1.max(seg_id);

            let len = seq.len();
            seg_lens.push(len);
        }

        assert!(
        seg_id_range.1 - seg_id_range.0 == seg_lens.len() - 1,
        "GFA segments must be tightly packed: min ID {}, max ID {}, node count {}, was {}",
        seg_id_range.0, seg_id_range.1, seg_lens.len(),
        seg_id_range.1 - seg_id_range.0,
        );

        let gfa = std::fs::File::open(&gfa_path)?;
        let mut gfa_reader = BufReader::new(gfa);

        let mut path_names = BTreeMap::default();

        let mut path_steps: Vec<Vec<PathStep>> = Vec::new();
        let mut path_step_offsets: Vec<RoaringBitmap> = Vec::new();
        // let mut path_pos: Vec<Vec<usize>> = Vec::new();

        loop {
            line_buf.clear();

            let len = gfa_reader.read_until(b'\n', &mut line_buf)?;
            if len == 0 {
                break;
            }

            let line = &line_buf[..len];
            if !matches!(line.first(), Some(b'P')) {
                continue;
            }

            let mut fields = line.split(|&c| c == b'\t');

            let Some((name, steps)) = fields.next().and_then(|_type| {
                let name = fields.next()?;
                let steps = fields.next()?;
                Some((name, steps))
            }) else {
                continue;
            };

            let name = std::str::from_utf8(name)?;
            path_names.insert(name.to_string(), path_steps.len());

            let mut pos = 0;

            let mut parsed_steps = Vec::new();

            let mut offsets = RoaringBitmap::new();

            let steps = steps.split(|&c| c == b',');

            for step in steps {
                let (seg, orient) = step.split_at(step.len() - 1);
                let seg_id = btoi::btou::<usize>(seg)?;
                let seg_ix = seg_id - seg_id_range.0;
                let len = seg_lens[seg_ix];

                let is_rev = orient == b"+";

                let step = PathStep {
                    node: seg_ix as u32,
                    reverse: is_rev,
                };
                parsed_steps.push(step);
                offsets.push(pos as u32);

                pos += len;
            }

            path_steps.push(parsed_steps);
            path_step_offsets.push(offsets);
        }

        Ok(Self {
            path_names,
            path_steps,
            path_step_offsets,

            segment_id_range: seg_id_range,
            segment_lens: seg_lens,
        })
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    let path_index = PathIndex::from_gfa(&args.gfa)?;

    if let Some(bam_path) = args.bam {
        return bam_injection(path_index, bam_path);
    } else if let Some(paf_path) = args.paf {
        return paf_injection(path_index, paf_path);
    } else if let Some(gbam_path) = args.gbam {
        return gbam_injection(path_index, gbam_path);
    } else if let Some((path, start, end)) = args.range {
        return path_range_cmd(path_index, path, start, end);
    }

    Ok(())
}

struct SamFlagInfo {
    is_paired: bool,
    is_first: bool,
    is_second: bool,
}

impl SamFlagInfo {
    fn from_flag(flag: u16) -> Self {
        SamFlagInfo {
            is_paired: (flag & 0x1) != 0,
            is_first: (flag & 0x40) != 0,
            is_second: (flag & 0x80) != 0,
        }
    }
}

fn process_query_name(record: &noodles::sam::alignment::record::Record) -> Option<String> {
    let read_name = record.read_name()?.to_string();
    let flags = SamFlagInfo::from_flag(record.flags().bits());
    
    // Return unmodified name if not paired
    if !flags.is_paired {
        return Some(read_name);
    }
    
    // Add appropriate suffix based on first/second read
    if flags.is_first {
        Some(format!("{}_R1", read_name))
    } else if flags.is_second {
        Some(format!("{}_R2", read_name))
    } else {
        Some(read_name)
    }
}

fn bam_injection(path_index: PathIndex, bam_path: PathBuf) -> Result<()> {
    use noodles::bam;

    let mut bam = std::fs::File::open(&bam_path).map(bam::Reader::new)?;

    let header = {
        use noodles::sam;
        // the noodles parse() impl demands that the @HD lines go first,
        // but that's clearly not a guarantee i can enforce
        let raw = bam.read_header()?;
        let mut header = sam::Header::builder();

        for line in raw.lines() {
            use noodles::sam::header::Record as HRecord;
            if let Ok(record) = line.parse::<HRecord>() {
                header = match record {
                    HRecord::Header(hd) => header.set_header(hd),
                    HRecord::ReferenceSequence(sq) => {
                        header.add_reference_sequence(sq)
                    }
                    HRecord::ReadGroup(rg) => header.add_read_group(rg),
                    HRecord::Program(pg) => header.add_program(pg),
                    HRecord::Comment(co) => header.add_comment(co),
                };
            }
        }

        header.build()
    };

    let _ref_seqs = bam.read_reference_sequences()?;

    // for (key, val) in ref_seqs {
    //     let len = val.length();
    //     eprintln!("{key}\t{len}");
    // }

    let mut stdout = std::io::stdout().lock();

    for rec in bam.records() {
        let record = rec?;

        // if record.flags().is_reverse_complemented() {
        //     continue;
        // }

        // Process read name with flags
        let Some(read_name) = process_query_name(&record) else {
            continue; // Skip record if name processing fails
        };

        // let name = read_name.to_string();
        // dbg!(&name);
        // let name = std::str::from_utf8(read_name.to_)?;

        // convenient for debugging
        // let name: &str = read_name.as_ref();
        // if name != "A00404:156:HV37TDSXX:4:1213:6370:23359" {
        //     continue;
        // }

        let Some(ref_name) = record
                    .reference_sequence(&header)
                    .and_then(|s| s.ok().map(|s| s.name()))
        else {
            continue;
        };

        let Some(path_id) = path_index.path_names.get(ref_name.as_str()).copied() else {
            continue;
        };

        // skip the record if there is no alignment information
        if record.alignment_start().is_none() || record.alignment_end().is_none() {
            continue;
        }

        let start = record.alignment_start().unwrap();
        let end = record.alignment_end().unwrap();
        let al_len = record.alignment_span();
        //assert!(end - start == al_len);

        let start_pos = (start.get()-1) as u32;
        let start_rank = path_index.path_step_offsets[path_id].rank(start_pos);
        //eprintln!("start_rank = {}", start_rank);
        let mut step_offset = start_pos
            - path_index.path_step_offsets[path_id]
                .select((start_rank - 1) as u32)
                .unwrap();

        let pos_range = ((start.get()-1) as u32)..((end.get()-1) as u32);
        if let Some(steps) =
            path_index.path_step_range_iter(ref_name.as_str(), pos_range)
        {
            let mut path_str = String::new();

            let mut steps = steps.collect::<Vec<_>>();

            if record.flags().is_reverse_complemented() {
                steps.reverse();
            }

            let mut path_len: usize = 0;

            for (_step_ix, step) in steps {
                // path length is given by the length of nodes in the graph
                path_len += path_index.segment_lens[(step.node) as usize];
                use std::fmt::Write;
                // eprintln!("step_ix: {step_ix}");

                // let reverse = step.reverse;
                let forward = step.reverse ^ record.flags().is_reverse_complemented();
                if forward {
                    write!(&mut path_str, ">")?;
                } else {
                    write!(&mut path_str, "<")?;
                }
                write!(
                    &mut path_str,
                    "{}",
                    step.node + path_index.segment_id_range.0 as u32
                )?;
            }

            if record.flags().is_reverse_complemented() {
                // start node offset changes
                //println!("is rev {} {} {}", path_len, step_offset, record.cigar().alignment_span());
                let last_bit = path_len as u32 - (step_offset as u32 + record.cigar().alignment_span() as u32 - 1);
                step_offset = last_bit;
            }

            // query name
            write!(stdout, "{}\t", read_name)?;

            // query len
            let query_len = record.cigar().read_length();
            write!(stdout, "{query_len}\t")?;

            // query start (0-based, closed)
            let query_start = 0;
            write!(stdout, "{query_start}\t")?;

            // query end (0-based, open)
            write!(stdout, "{}\t", query_start + query_len)?;

            // strand
            // if record.flags().is_reverse_complemented() {
            // print!("-\t");
            // } else {
            write!(stdout, "+\t")?;
            // }

            // path
            write!(stdout, "{path_str}\t")?;
            // path length
            write!(stdout, "{path_len}\t")?;
            // start on path
            let path_start = step_offset as usize;
            write!(stdout, "{path_start}\t")?;
            // end on path
            let path_end = path_start + al_len;
            write!(stdout, "{path_end}\t")?;
            // number of matches
            {
                use noodles::sam::record::cigar::{op::Kind, Op};

                fn match_len(op: &Op) -> usize {
                    match op.kind() {
                        Kind::Match
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch => op.len(),
                        _ => 0,
                    }
                }
                let matches =
                    record.cigar().iter().map(match_len).sum::<usize>();
                write!(stdout, "{matches}\t")?;
            }
          // alignment block length
            write!(stdout, "{al_len}\t")?;
            // mapping quality
            {
                let score =
                    record.mapping_quality().map(|q| q.get()).unwrap_or(255u8);
                write!(stdout, "{score}\t")?;
            }

            // cigar
            write!(stdout, "cg:Z:{}", record.cigar())?;

            writeln!(stdout)?;
        } else {
        }
    }

    std::io::stdout().flush()?;

    Ok(())
}

fn gbam_injection(path_index: PathIndex, gbam_path: PathBuf) -> Result<()> {
    let file = File::open(gbam_path.clone()).unwrap();
    let mut template = ParsingTemplate::new();
    // Only fetch fields which are needed.
    template.set(&Fields::Flags, true);
    template.set(&Fields::RefID, true);
    template.set(&Fields::ReadName, true);
    template.set(&Fields::RawCigar, true);
    template.set(&Fields::Mapq, true);
    template.set(&Fields::Pos, true);

    let mut reader = Reader::new(file, template).unwrap();

    let ref_seqs = reader.file_meta.get_ref_seqs().clone();
    // dbg!(&ref_seqs);
    // return Ok(());

    let mut records_it = Records::new(&mut reader);

    let mut stdout = std::io::stdout().lock();

    while let Some(rec) = records_it.next_rec() {
        if rec.is_unmapped() || rec.refid.unwrap() < 0 {
            continue; // Unmapped read
        }
        let ref_name = &ref_seqs[rec.refid.unwrap() as usize].0;
        let Some(path_id) = path_index.path_names.get(ref_name.as_str()).copied() else {
            continue;
        };

        // skip the record if there is no alignment information
        if rec.alignment_start().is_none() || rec.alignment_end().is_none() {
            continue;
        }

        let start = rec.alignment_start().unwrap();
        let end = rec.alignment_end().unwrap();
        let al_len = rec.alignment_span() as usize;
        //assert!(end - start == al_len);

        let start_pos = start as u32;
        let start_rank = path_index.path_step_offsets[path_id].rank(start_pos);
        //eprintln!("start_rank = {}", start_rank);
        let mut step_offset = start_pos
            - path_index.path_step_offsets[path_id]
            .select((start_rank - 1) as u32)
            .unwrap();

        let pos_range = ((start) as u32)..((end) as u32);
        if let Some(steps) =
            path_index.path_step_range_iter(ref_name.as_str(), pos_range)
        {
            let mut path_str = String::new();

            let mut steps = steps.collect::<Vec<_>>();

            if rec.is_reverse_complemented() {
                steps.reverse();
            }

            let mut path_len: usize = 0;

            for (_step_ix, step) in steps {
                // path length is given by the length of nodes in the graph
                path_len += path_index.segment_lens[(step.node) as usize];
                use std::fmt::Write;
                // eprintln!("step_ix: {step_ix}");

                // let reverse = step.reverse;
                let forward = step.reverse ^ rec.is_reverse_complemented();
                if forward {
                    write!(&mut path_str, ">")?;
                } else {
                    write!(&mut path_str, "<")?;
                }
                write!(
                    &mut path_str,
                    "{}",
                    step.node + path_index.segment_id_range.0 as u32
                )?;
            }

            if rec.is_reverse_complemented() {
                // start node offset changes
                //println!("is rev {} {} {}", path_len, step_offset, record.cigar().alignment_span());
                let last_bit = path_len as u32 - (step_offset as u32 + rec.alignment_span() as u32 - 1);
                step_offset = last_bit;
            }

            // query name
            // let read_name =  String::from_utf8(rec.read_name.clone().unwrap()).unwrap();
            let read_name = unsafe {std::str::from_utf8_unchecked(&rec.read_name.as_ref().unwrap())};

            write!(stdout, "{}\t", read_name.trim_end_matches('\0'))?;

            // query len
            let query_len = rec.cigar.as_ref().unwrap().read_length();
            write!(stdout, "{query_len}\t")?;

            //todo to check!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // query start (0-based, closed)
            let query_start = 0;
            write!(stdout, "{query_start}\t")?;

            // query end (0-based, open)
            write!(stdout, "{}\t", query_start + query_len)?;

            //todo to check!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // strand
            // if rec.is_reverse_complemented() {
            // print!("-\t");
            // } else {
            write!(stdout, "+\t")?;
            // }

            // path
            write!(stdout, "{path_str}\t")?;
            // path length
            write!(stdout, "{path_len}\t")?;
            // start on path
            let path_start = step_offset as usize;
            write!(stdout, "{path_start}\t")?;
            // end on path
            let path_end = path_start + al_len;
            write!(stdout, "{path_end}\t")?;
            // number of matches
            {
                fn match_len(op: &gbam_tools::query::cigar::Op) -> usize {
                    match op.op_type() {
                        'M' | '=' | 'X' => op.length() as usize,
                        _ => 0,
                    }
                }
                let matches =
                    rec.cigar.as_ref().unwrap().ops().map(match_len).sum::<usize>();
                //println!("\nmatches: {}", matches);
                write!(stdout, "{matches}\t")?;
            }
            // alignment block length
            write!(stdout, "{al_len}\t")?;
            // mapping quality
            {
                let score =
                    rec.mapq.map(|q| q).unwrap_or(255u8);
                write!(stdout, "{score}\t")?;
            }

            // cigar
            write!(stdout, "cg:Z:{}",rec.cigar.as_ref().unwrap())?;

            writeln!(stdout)?;
        }
    }

    std::io::stdout().flush()?;

    Ok(())
}

fn paf_injection(path_index: PathIndex, paf_path: PathBuf) -> Result<()> {
    let file = std::fs::File::open(paf_path)?;
    let reader = BufReader::new(file);
    let mut stdout = std::io::stdout().lock();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            continue; // Skip invalid lines
        }

        let query_name = fields[0];
        let query_len: usize = fields[1].parse()?;
        let query_start: usize = fields[2].parse()?;
        let query_end: usize = fields[3].parse()?;
        let strand = fields[4];
        let ref_name = fields[5];
        let _ref_len: usize = fields[6].parse()?;
        let ref_start: usize = fields[7].parse()?;
        let _ref_end: usize = fields[8].parse()?;
        let mapping_quality: u8 = fields[11].parse()?;

        // Find the cg:Z: tag for CIGAR string
        let cigar_str = fields.iter().find(|&&f| f.starts_with("cg:Z:"))
            .map(|&f| &f[5..])
            .unwrap_or("");
        let (alignment_span, matches) = calculate_alignment_stats(cigar_str);


        if let Some(path_id) = path_index.path_names.get(ref_name).copied() {
            let pos_range = (ref_start as u32)..((ref_start + alignment_span - 1) as u32);
        
            if let Some(steps) = path_index.path_step_range_iter(ref_name, pos_range) {
                let mut path_str = String::new();
                let mut path_len: usize = 0;

                let is_reverse_complemented = strand == "-";

                let mut steps = steps.collect::<Vec<_>>();

                if is_reverse_complemented {
                    steps.reverse();
                }

                for (_step_ix, step) in steps {
                    path_len += path_index.segment_lens[(step.node) as usize];
                    let forward = step.reverse ^ is_reverse_complemented;
                    use std::fmt::Write;
                    write!(&mut path_str, "{}{}", if forward { ">" } else { "<" }, step.node + path_index.segment_id_range.0 as u32)?;
                }

                // Calculate step_offset
                let start_pos = ref_start as u32;
                let start_rank = path_index.path_step_offsets[path_id].rank(start_pos);
                let mut step_offset = start_pos
                    - path_index.path_step_offsets[path_id]
                        .select((start_rank - 1) as u32)
                        .unwrap();

                if is_reverse_complemented {
                    let last_bit = path_len as u32 - (step_offset as u32 + alignment_span as u32 - 1);
                    step_offset = last_bit;
                }

                // Calculate path_start and path_end
                let path_start = step_offset as usize;
                let path_end = path_start + alignment_span;

                // Output in the same format as the BAM processing
                writeln!(stdout, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}",
                    query_name, query_len, query_start, query_end, "+",
                    path_str, path_len, path_start, path_end, matches, alignment_span, mapping_quality, cigar_str)?;
            }
        }
    }

    Ok(())
}

/// Calculate alignment span and number of matches from a CIGAR string
/// Returns (total_span, num_matches)
fn calculate_alignment_stats(cigar: &str) -> (usize, usize) {
    let mut total_span = 0;
    let mut num_matches = 0;
    let mut num_buffer = String::with_capacity(4);

    // Parse each character in the CIGAR string
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_buffer.push(ch);
        } else if !num_buffer.is_empty() {
            let count: usize = num_buffer.parse().unwrap();
            // Update span for alignment operations
            match ch {
                'M' | 'D' | 'N' | '=' | 'X' => total_span += count,
                _ => {}
            }
            
            // Count matches separately
            if ch == '=' || ch == 'M' {
                num_matches += count;
            }
            num_buffer.clear();
        }
    }

    (total_span, num_matches)
}

fn path_range_cmd(
    path_index: PathIndex,
    path_name: String,
    start: usize,
    end: usize,
) -> Result<()> {
    let path = path_index
        .path_names
        .get(&path_name)
        .expect("Path not found");

    let offsets = path_index.path_step_offsets.get(*path).unwrap();
    let steps = path_index.path_steps.get(*path).unwrap();

    let start_rank = offsets.rank(start as u32);
    let end_rank = offsets.rank(end as u32);

    let cardinality = offsets.range_cardinality((start as u32)..(end as u32));

    println!("start_rank: {start_rank}");
    println!("end_rank: {end_rank}");
    println!("cardinality: {cardinality}");

    println!("------");
    let skip = (start_rank as usize).checked_sub(1).unwrap_or_default();
    let take = end_rank as usize - skip;


    println!("step_ix\tnode\tpos");
    // let skip = 0;
    // let take = 20;
    for (step_ix, (pos, step)) in
        offsets.iter().zip(steps).enumerate().skip(skip).take(take)
    {
        let node = step.node + path_index.segment_id_range.0 as u32;
        println!("{step_ix}\t{node}\t{pos}");
    }

    println!("------------");

    // let start_pos = start.get() as u32;
    // let start_rank = path_index.path_step_offsets[path_id].rank(start_pos);
    // let step_offset = start_pos
    //     - path_index.path_step_offsets[path_id]
    //         .select((start_rank - 1) as u32)
    //         .unwrap();

    // let pos_range = (start.get() as u32)..(1 + end.get() as u32);
    let pos_range = (start as u32)..(end as u32);
    if let Some(steps) = path_index.path_step_range_iter(&path_name, pos_range)
    {
        for (step_ix, step) in steps {
            let pos = offsets.select(step_ix as u32).unwrap();
            //
            let node = step.node + path_index.segment_id_range.0 as u32;
            println!("{step_ix}\t{node}\t{pos}");
        }
    }

    Ok(())
}
