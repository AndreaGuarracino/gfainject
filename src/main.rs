use anyhow::Result;
use roaring::RoaringBitmap;
use std::collections::BTreeMap;
use std::io::prelude::*;
use std::{io::stdin, io::BufReader, path::PathBuf};

use gbam_tools::reader::parse_tmplt::ParsingTemplate;
use gbam_tools::reader::reader::Reader;
use gbam_tools::reader::records::Records;
use gbam_tools::Fields;
use std::fs::File;

use itertools::Itertools;

use clap::{ArgGroup, Parser};
use std::num::NonZeroUsize;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(group(
    ArgGroup::new("input_mode")
        .required(true)
        .args(["sam", "bam", "paf", "gbam", "range"]),
))]
struct Args {
    /// Path to input GFA file
    #[arg(long, required = true)]
    gfa: PathBuf,

    /// Path to input SAM file
    #[arg(short, long)]
    sam: Option<PathBuf>,

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

    /// Emit up to ALT_HITS alternative alignments (from XA tag)
    #[arg(long)]
    alt_hits: Option<NonZeroUsize>,
}

/// Parse a range string in the format "path_name:start-end"
/// where path_name can contain ':' and '-' characters
fn parse_range(s: &str) -> Result<(String, usize, usize), String> {
    // Find the last colon in the string
    let last_colon_pos = s
        .rfind(':')
        .ok_or("Range must include ':' separator".to_string())?;

    // Split into path_name and range parts
    let path_name = s[..last_colon_pos].to_string();
    let range_str = &s[last_colon_pos + 1..];

    // Split the range part by the last hyphen
    let range_parts: Vec<&str> = range_str.split('-').collect();
    if range_parts.len() != 2 {
        return Err("Range must be in format start-end".to_string());
    }

    // Parse start and end positions
    let start = range_parts[0]
        .parse::<usize>()
        .map_err(|_| "Invalid start position".to_string())?;
    let end = range_parts[1]
        .parse::<usize>()
        .map_err(|_| "Invalid end position".to_string())?;

    if end <= start {
        return Err(
            "End position must be greater than start position".to_string()
        );
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
    // path_id: usize,
    // pos_range: std::ops::Range<u32>,
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
            // path_id,
            // pos_range,
            steps,
            // first_step_start_pos,
            // last_step_end_pos,
        })
    }

    fn from_gfa(gfa_path: impl AsRef<std::path::Path>) -> Result<Self> {
        // First pass: Read segments to get lengths and ID range
        let gfa = std::fs::File::open(&gfa_path)?;
        let mut gfa_reader = BufReader::new(gfa);
        let mut line_buf = Vec::new();
        let mut seg_lens = Vec::new();

        let mut seg_id_range = (usize::MAX, 0usize);
        // dbg!();

        loop {
            line_buf.clear();
            let len = gfa_reader.read_until(b'\n', &mut line_buf)?;
            if len == 0 {
                break; // End of file
            }

            let line = &line_buf[..len];
            let line_str = std::str::from_utf8(line)?;
            // println!("{line_str}");

            // Skip non-segment lines
            if !matches!(line.first(), Some(b'S')) {
                continue;
            }

            // Parse segment line format: S<tab>id<tab>sequence
            let mut fields = line_str.split('\t');
            let Some((name, seq)) = fields.next().and_then(|_type| {
                let name = fields.next()?.trim();
                let seq = fields.next()?.trim();
                Some((name, seq))
            }) else {
                continue;
            };

            // Track segment ID range and store sequence length
            let seg_id = name.parse::<usize>()?;
            seg_id_range.0 = seg_id_range.0.min(seg_id);
            seg_id_range.1 = seg_id_range.1.max(seg_id);
            let len = seq.len();
            seg_lens.push(len);
        }

        // Verify segments are tightly packed (no gaps in ID numbering)
        assert!(
            seg_id_range.1 - seg_id_range.0 == seg_lens.len() - 1,
            "GFA segments must be tightly packed: min ID {}, max ID {}, node count {}, was {}",
            seg_id_range.0, seg_id_range.1, seg_lens.len(),
            seg_id_range.1 - seg_id_range.0,
        );

        // Second pass: Read paths to build path information
        let gfa = std::fs::File::open(&gfa_path)?;
        let mut gfa_reader = BufReader::new(gfa);
        let mut path_names = BTreeMap::default();
        let mut path_steps: Vec<Vec<PathStep>> = Vec::new();
        let mut path_step_offsets: Vec<RoaringBitmap> = Vec::new();
        // let mut path_pos: Vec<Vec<usize>> = Vec::new();

        // Process each line looking for Path ('P') records
        loop {
            line_buf.clear();
            let len = gfa_reader.read_until(b'\n', &mut line_buf)?;
            if len == 0 {
                break; // End of file
            }

            let line = &line_buf[..len];
            if !matches!(line.first(), Some(b'P')) {
                continue;
            }

            // Parse path line format: P<tab>name<tab>segments
            let mut fields = line.split(|&c| c == b'\t');
            let Some((name, steps)) = fields.next().and_then(|_type| {
                let name = fields.next()?;
                let steps = fields.next()?;
                Some((name, steps))
            }) else {
                continue;
            };

            // Store path name and its index
            let name = std::str::from_utf8(name)?;
            path_names.insert(name.to_string(), path_steps.len());

            let mut pos = 0; // Track position along the path
            let mut parsed_steps = Vec::new();
            let mut offsets = RoaringBitmap::new();

            // Process comma-separated path steps
            let steps = steps.split(|&c| c == b',');
            for step in steps {
                let (seg, orient) = step.split_at(step.len() - 1);
                let seg_id = btoi::btou::<usize>(seg)?;
                let seg_ix = seg_id - seg_id_range.0;
                let len = seg_lens[seg_ix];

                let is_rev = orient == b"-";

                // Create and store path step
                let step = PathStep {
                    node: seg_ix as u32,
                    reverse: is_rev,
                };
                parsed_steps.push(step);
                let _ = offsets.try_push(pos as u32);

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

    if let Some(sam_path) = args.sam {
        return sam_injection(path_index, sam_path, args.alt_hits);
    } else if let Some(bam_path) = args.bam {
        return bam_injection(path_index, bam_path, args.alt_hits);
    } else if let Some(paf_path) = args.paf {
        return paf_injection(path_index, paf_path, args.alt_hits);
    } else if let Some(gbam_path) = args.gbam {
        return gbam_injection(path_index, gbam_path, args.alt_hits);
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

#[derive(Debug)]
struct AlternativeHit {
    chr: String,
    strand: bool, // true for forward (+), false for reverse (-)
    pos: u32,
    cigar: String,
    nm: u32,
}

impl AlternativeHit {
    fn from_xa_str(xa_str: &str) -> Option<Self> {
        let parts: Vec<&str> = xa_str.split(',').collect();
        if parts.len() != 4 {
            return None;
        }

        let chr = parts[0].to_string();
        let pos_str = parts[1];
        let (strand, pos) = if let Some(stripped) = pos_str.strip_prefix('-') {
            (false, stripped.parse::<u32>().ok()?)
        } else if let Some(stripped) = pos_str.strip_prefix('+') {
            (true, stripped.parse::<u32>().ok()?)
        } else {
            return None;
        };
        let cigar = parts[2].to_string();
        let nm = parts[3].parse::<u32>().ok()?;

        Some(AlternativeHit {
            chr,
            strand,
            pos,
            cigar,
            nm,
        })
    }
}

fn process_query_name(query_name: &str, flags: SamFlagInfo) -> String {
    // Return unmodified name if not paired
    if !flags.is_paired {
        return query_name.to_string();
    }

    // Add appropriate suffix based on first/second read
    if flags.is_first {
        format!("{}/1", query_name)
    } else if flags.is_second {
        format!("{}/2", query_name)
    } else {
        query_name.to_string()
    }
}

// Function to write GAF format output
fn write_gaf_record<W: std::io::Write>(
    writer: &mut W,
    read_name: &str,
    query_len: usize,
    query_start: usize,
    query_end: usize,
    path_str: &str,
    path_len: usize,
    path_start: usize,
    path_end: usize,
    matches: usize,
    alignment_span: usize,
    mapping_quality: u8,
    cigar: &str,
) -> Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t+\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}",
        read_name,
        query_len,
        query_start,
        query_end,
        path_str,
        path_len,
        path_start,
        path_end,
        matches,
        alignment_span,
        mapping_quality,
        cigar
    )?;
    Ok(())
}

fn calculate_query_coords_and_align_stats(
    cigar: &str,
) -> (usize, usize, usize, usize, usize) {
    let mut query_start = 0;
    let mut query_len = 0;
    let mut query_consumed = 0; // Track query bases consumed by M/=/X/I operations
    let mut total_span = 0; // Total span of alignment operations that consume the reference sequence (M/D/N/=/X operations)
    let mut num_matches = 0;
    let mut num_buffer = String::with_capacity(4);
    let mut found_first_match = false;

    // Parse each character in the CIGAR string
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num_buffer.push(ch);
        } else if !num_buffer.is_empty() {
            let count: usize = num_buffer.parse().unwrap();

            // 'S', 'H', 'P', 'I' do not consume the reference sequence
            match ch {
                'S' | 'H' => {
                    if !found_first_match {
                        query_start += count;
                    }

                    query_len += count;
                }
                'M' | '=' | 'X' => {
                    found_first_match = true;

                    query_len += count;
                    query_consumed += count;
                    total_span += count;
                    num_matches += count;
                }
                'D' | 'N' => total_span += count,
                'I' => {
                    found_first_match = true;

                    query_len += count;
                    query_consumed += count;
                }
                _ => {} // 'P' and others
            }
            num_buffer.clear();
        }
    }

    (
        query_len,
        query_start,
        query_start + query_consumed,
        total_span,
        num_matches,
    )
}

#[derive(Debug)]
struct AlignmentInfo {
    ref_name: String,
    read_name: String,
    ref_start_pos: u32,
    is_reverse: bool,
    cigar_str: String,
    mapping_quality: u8,
}

fn process_alignment(
    alignment: AlignmentInfo,
    path_index: &PathIndex,
    stdout: &mut impl std::io::Write,
) -> Result<()> {
    if alignment.cigar_str.is_empty() || alignment.cigar_str == "*" {
        return Ok(());
    }

    let Some(path_id) = path_index.path_names.get(&alignment.ref_name).copied()
    else {
        return Ok(());
    };

    let (query_len, query_start, query_end, alignment_span, num_matches) =
        calculate_query_coords_and_align_stats(&alignment.cigar_str);
    let alignment_end_pos =
        alignment.ref_start_pos + (alignment_span - 1) as u32;
    let pos_range = alignment.ref_start_pos..alignment_end_pos;

    if let Some(steps) =
        path_index.path_step_range_iter(&alignment.ref_name, pos_range)
    {
        use std::fmt::Write;

        let mut path_str = String::new();
        let mut steps = steps.collect::<Vec<_>>();
        let mut path_len: usize = 0;

        if alignment.is_reverse {
            steps.reverse();
        }

        for (_step_ix, step) in steps {
            // path length is given by the length of nodes in the graph
            path_len += path_index.segment_lens[step.node as usize];
            // eprintln!("step_ix: {step_ix}");

            // let forward = step.reverse ^ alignment.is_reverse;
            let forward = step.reverse == alignment.is_reverse;

            write!(
                &mut path_str,
                "{}{}",
                if forward { ">" } else { "<" },
                step.node + path_index.segment_id_range.0 as u32
            )?;
        }

        // Find starting position within first node
        let start_rank =
            path_index.path_step_offsets[path_id].rank(alignment.ref_start_pos);
        //eprintln!("start_rank = {}", start_rank);
        let mut step_offset = alignment.ref_start_pos
            - path_index.path_step_offsets[path_id]
                .select((start_rank - 1) as u32)
                .unwrap();
        // Adjust offset for reverse alignments
        if alignment.is_reverse {
            // start node offset changes
            //let last_bit = path_len as u32 - (step_offset as u32 + record.cigar().alignment_span() as u32 - 1);
            let last_bit =
                path_len as u32 - (step_offset + alignment_span as u32);
            step_offset = last_bit;
        }

        write_gaf_record(
            stdout,
            &alignment.read_name,
            query_len,
            query_start,
            query_end,
            &path_str,
            path_len,
            step_offset as usize,
            step_offset as usize + alignment_span,
            num_matches,
            alignment_span,
            alignment.mapping_quality,
            &alignment.cigar_str,
        )?;
    }

    Ok(())
}

fn sam_injection(
    path_index: PathIndex,
    sam_path: PathBuf,
    alt_hits: Option<NonZeroUsize>,
) -> Result<()> {
    use std::io::{BufRead, BufReader};

    // Support stdin when path is "-"
    let reader: Box<dyn BufRead> = if sam_path.as_os_str() == "-" {
        Box::new(BufReader::new(stdin()))
    } else {
        Box::new(BufReader::new(std::fs::File::open(sam_path)?))
    };
    let mut stdout = std::io::stdout().lock();

    for line in reader.lines() {
        let line = line?;

        // Skip header lines
        if line.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        // Parse SAM fields
        let qname = fields[0];
        let flag = fields[1].parse::<u16>()?;
        let rname = fields[2];
        let pos = fields[3].parse::<u32>()? - 1; // Convert to 0-based
        let mapq = fields[4].parse::<u8>()?;
        let cigar_str = fields[5];

        // Skip unmapped reads
        if flag & 0x4 != 0 || rname == "*" || cigar_str == "*" {
            continue;
        }

        // Skip secondary alignments
        if flag & 0x100 != 0 {
            continue;
        }

        // Process read name with flags
        let flags = SamFlagInfo::from_flag(flag);
        let read_name = process_query_name(qname, flags);

        // Process primary alignment
        let primary_alignment = AlignmentInfo {
            ref_name: rname.to_string(),
            read_name: read_name.clone(),
            ref_start_pos: pos,
            is_reverse: flag & 0x10 != 0,
            cigar_str: cigar_str.to_string(),
            mapping_quality: mapq,
        };
        process_alignment(primary_alignment, &path_index, &mut stdout)?;

        // Process alternative alignments if requested
        let Some(max_alt) = alt_hits else {
            continue;
        };

        // Find XA and NM tags
        let mut xa_str = None;
        let mut primary_nm = None;

        for field in fields.iter().skip(11) {
            if let Some(stripped) = field.strip_prefix("XA:Z:") {
                xa_str = Some(stripped);
            } else if let Some(stripped) = field.strip_prefix("NM:i:") {
                primary_nm = stripped.parse::<u32>().ok();
            }
        }

        let (Some(xa), Some(nm)) = (xa_str, primary_nm) else {
            continue;
        };

        // Process alternative alignments
        let alt_hits_vec: Vec<_> = xa
            .split(';')
            .filter(|hit| !hit.is_empty())
            .filter_map(AlternativeHit::from_xa_str)
            .sorted_by_key(|hit| hit.nm)
            .take_while(|hit| hit.nm <= nm)
            .take(max_alt.into())
            .collect();

        for alt_hit in alt_hits_vec {
            let alt_alignment = AlignmentInfo {
                ref_name: alt_hit.chr,
                read_name: read_name.clone(),
                ref_start_pos: alt_hit.pos - 1,
                is_reverse: !alt_hit.strand,
                cigar_str: alt_hit.cigar,
                mapping_quality: 255,
            };
            process_alignment(alt_alignment, &path_index, &mut stdout)?;
        }
    }

    Ok(())
}

fn bam_injection(
    path_index: PathIndex,
    bam_path: PathBuf,
    alt_hits: Option<NonZeroUsize>,
) -> Result<()> {
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

    // IMPORTANT: this instruction has to be execute to avoid 'Error: invalid read name terminator: expected 0x00, got 0x65'
    let _ref_seqs = bam.read_reference_sequences()?;
    // for (key, val) in _ref_seqs {
    //     let len = val.length();
    //     eprintln!("{key}\t{len}");
    // }

    let mut stdout = std::io::stdout().lock();

    for rec in bam.records() {
        let record = rec?;

        // skip the record if there is no alignment information
        let (Some(start), Some(_end)) =
            (record.alignment_start(), record.alignment_end())
        else {
            continue;
        };

        let Some(ref_name) = record
            .reference_sequence(&header)
            .and_then(|s| s.ok().map(|s| s.name()))
        else {
            continue;
        };

        // Process read name with flags
        let Some(read_name) = record.read_name() else {
            continue;
        };
        let flags = SamFlagInfo::from_flag(record.flags().bits());
        let read_name = process_query_name(read_name.as_ref(), flags);

        // let name = read_name.to_string();
        // dbg!(&name);

        // convenient for debugging
        // let name: &str = read_name.as_ref();
        // if name != "A00404:156:HV37TDSXX:4:1213:6370:23359" {
        //     continue;
        // }

        let primary_alignment = AlignmentInfo {
            ref_name: ref_name.to_string(),
            read_name: read_name.clone(),
            //read_len: record.sequence().len(),
            ref_start_pos: (start.get() - 1) as u32,
            is_reverse: record.flags().is_reverse_complemented(),
            cigar_str: record.cigar().to_string(),
            mapping_quality: record
                .mapping_quality()
                .map(|q| q.get())
                .unwrap_or(255),
        };
        process_alignment(primary_alignment, &path_index, &mut stdout)?;

        // Process alternative alignments only if requested
        let Some(max_alt) = alt_hits else {
            continue;
        };

        let xa_str = {
            use noodles::sam::record::data::field::Tag;
            use std::str::FromStr;

            let xa_tag = Tag::from_str("XA").expect("Valid tag");
            if let Some(field) = record.data().get(xa_tag) {
                if let noodles::sam::record::data::field::Value::String(
                    xa_str,
                ) = field.value()
                {
                    Some(xa_str)
                } else {
                    None
                }
            } else {
                None
            }
        };

        let Some(xa_str) = xa_str else {
            continue;
        };

        // Get NM tag from primary alignment if alt_hits enabled
        let primary_nm = {
            use noodles::sam::record::data::field::Tag;
            use std::str::FromStr;

            let nm_tag = Tag::from_str("NM").expect("Valid tag");
            if let Some(field) = record.data().get(nm_tag) {
                match field.value() {
                    noodles::sam::record::data::field::Value::UInt8(nm) => {
                        Some(*nm as u32)
                    }
                    noodles::sam::record::data::field::Value::UInt16(nm) => {
                        Some(*nm as u32)
                    }
                    noodles::sam::record::data::field::Value::UInt32(nm) => {
                        Some(*nm)
                    }
                    _ => None,
                }
            } else {
                None
            }
        };

        let Some(primary_nm) = primary_nm else {
            continue;
        };

        // Collect all alternative hits and sort by NM value
        let alt_hits_vec: Vec<_> = xa_str
            .split(';')
            .filter(|hit| !hit.is_empty())
            .filter_map(AlternativeHit::from_xa_str)
            .sorted_by_key(|hit| hit.nm)
            .take_while(|hit| hit.nm <= primary_nm)
            .take(max_alt.into())
            .collect();

        // Take up to max_alt hits with NM <= primary_nm
        for alt_hit in alt_hits_vec.into_iter() {
            let alt_alignment = AlignmentInfo {
                ref_name: alt_hit.chr,
                read_name: read_name.clone(),
                ref_start_pos: alt_hit.pos - 1,
                is_reverse: !alt_hit.strand,
                cigar_str: alt_hit.cigar,
                mapping_quality: 255,
            };
            process_alignment(alt_alignment, &path_index, &mut stdout)?;
        }
    }

    std::io::stdout().flush()?;

    Ok(())
}

// Add function to parse XA tag in GBAM
fn parse_xa_tag(tags: &[u8]) -> Option<String> {
    let mut i = 0;
    while i < tags.len() - 2 {
        // Look for "XA" followed by tag type identifier
        if tags[i] == b'X' && tags[i + 1] == b'A' && tags[i + 2] == b'Z' {
            // Skip tag name and type
            i += 3;
            let mut result = Vec::new();

            // Read until null terminator or end of tags
            while i < tags.len() && tags[i] != 0 {
                result.push(tags[i]);
                i += 1;
            }

            if !result.is_empty() {
                return String::from_utf8(result).ok();
            }
        }
        i += 1;
    }
    None
}

// Add function to parse NM tag in GBAM
fn parse_nm_tag(tags: &[u8]) -> Option<u32> {
    let mut i = 0;
    while i < tags.len() - 2 {
        // Look for "NM"
        if tags[i] == b'N' && tags[i + 1] == b'M' {
            i += 2;
            // Skip type identifier if present
            if i < tags.len() && tags[i] == b'C' {
                i += 1;
            }
            // Get the value
            if i < tags.len() {
                return Some(tags[i] as u32);
            }
        }
        i += 1;
    }
    None
}

fn gbam_injection(
    path_index: PathIndex,
    gbam_path: PathBuf,
    alt_hits: Option<NonZeroUsize>,
) -> Result<()> {
    let file = File::open(gbam_path.clone()).unwrap();
    let mut template = ParsingTemplate::new();
    // Only fetch fields which are needed.
    template.set(&Fields::Flags, true);
    template.set(&Fields::RefID, true);
    template.set(&Fields::ReadName, true);
    template.set(&Fields::RawCigar, true);
    template.set(&Fields::Mapq, true);
    template.set(&Fields::Pos, true);
    if alt_hits.is_some() {
        template.set(&Fields::RawTags, true); // Need tags for XA and NM fields
    }

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

        // skip the record if there is no alignment information
        let (Some(start), Some(_end)) =
            (rec.alignment_start(), rec.alignment_end())
        else {
            continue;
        };

        let ref_name = &ref_seqs[rec.refid.unwrap() as usize].0;

        // Process read name with flags
        let read_name = unsafe {
            std::str::from_utf8_unchecked(rec.read_name.as_ref().unwrap())
        };
        let flags = SamFlagInfo::from_flag(rec.flag.unwrap());
        let read_name =
            process_query_name(read_name.trim_end_matches('\0'), flags);

        let primary_alignment = AlignmentInfo {
            ref_name: ref_name.to_string(),
            read_name: read_name.clone(),
            //read_len: rec.cigar.as_ref().unwrap().read_length() as usize,
            ref_start_pos: start, // already 0-based
            is_reverse: rec.is_reverse_complemented(),
            cigar_str: rec.cigar.as_ref().unwrap().to_string(),
            mapping_quality: rec.mapq.unwrap_or(255),
        };
        process_alignment(primary_alignment, &path_index, &mut stdout)?;

        // Process alternative alignments only if requested
        let Some(max_alt) = alt_hits else {
            continue;
        };

        // Get the XA string and primary NM value
        let (Some(xa_str), Some(primary_nm)) = (
            rec.tags.as_ref().and_then(|tags| parse_xa_tag(tags)),
            rec.tags.as_ref().and_then(|tags| parse_nm_tag(tags)),
        ) else {
            continue;
        };

        // Collect and sort alternative hits
        let alt_hits_vec: Vec<_> = xa_str
            .split(';')
            .filter(|hit| !hit.is_empty())
            .filter_map(AlternativeHit::from_xa_str)
            .sorted_by_key(|hit| hit.nm)
            .take_while(|hit| hit.nm <= primary_nm)
            .take(max_alt.into())
            .collect();

        // Take up to max_alt hits with NM <= primary_nm
        for alt_hit in alt_hits_vec.into_iter() {
            let alt_alignment = AlignmentInfo {
                ref_name: alt_hit.chr,
                read_name: read_name.clone(),
                //read_len: rec.cigar.as_ref().unwrap().read_length() as usize,
                ref_start_pos: alt_hit.pos - 1,
                is_reverse: !alt_hit.strand,
                cigar_str: alt_hit.cigar,
                mapping_quality: 255,
            };
            process_alignment(alt_alignment, &path_index, &mut stdout)?;
        }
    }

    std::io::stdout().flush()?;

    Ok(())
}

fn paf_injection(
    path_index: PathIndex,
    paf_path: PathBuf,
    alt_hits: Option<NonZeroUsize>,
) -> Result<()> {
    let file = std::fs::File::open(paf_path)?;
    let reader = BufReader::new(file);
    let mut stdout = std::io::stdout().lock();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            continue; // Skip invalid lines
        }

        // Find the cg:Z: tag for CIGAR string
        let cigar_str = fields
            .iter()
            .find(|&&f| f.starts_with("cg:Z:"))
            .map(|&f| &f[5..])
            .unwrap_or("");
        if cigar_str.is_empty() || cigar_str == "*" {
            continue;
        }

        // let query_len: usize = fields[1].parse()?;
        // let query_start: usize = fields[2].parse()?;
        // let query_end: usize = fields[3].parse()?;
        // let _ref_len: usize = fields[6].parse()?;
        // let _ref_end: usize = fields[8].parse()?;

        let alignment = AlignmentInfo {
            ref_name: fields[5].to_string(),
            read_name: fields[0].to_string(),
            //read_len: fields[1].parse()?,
            ref_start_pos: fields[7].parse::<u32>()?, // already 0-based
            is_reverse: fields[4] == "-",
            cigar_str: cigar_str.to_string(),
            mapping_quality: fields[11].parse()?,
        };
        process_alignment(alignment, &path_index, &mut stdout)?;

        // Process alternative alignments only if requested
        let Some(max_alt) = alt_hits else {
            continue;
        };

        // Find XA:Z: tag for alternative alignments
        let xa_str = fields
            .iter()
            .find(|&&f| f.starts_with("XA:Z:"))
            .map(|&f| &f[5..]);
        let Some(xa_str) = xa_str else {
            continue;
        };

        // Find NM:i: tag for primary alignment
        let primary_nm = fields
            .iter()
            .find(|&&f| f.starts_with("NM:i:"))
            .and_then(|&f| f[5..].parse::<u32>().ok());
        let Some(primary_nm) = primary_nm else {
            continue;
        };

        // Collect and sort alternative hits by NM value
        let alt_hits_vec: Vec<_> = xa_str
            .split(';')
            .filter(|hit| !hit.is_empty())
            .filter_map(AlternativeHit::from_xa_str)
            .sorted_by_key(|hit| hit.nm)
            .take_while(|hit| hit.nm <= primary_nm)
            .take(max_alt.into())
            .collect();

        // Process alternative alignments
        for alt_hit in alt_hits_vec.into_iter() {
            let alt_alignment = AlignmentInfo {
                ref_name: alt_hit.chr,
                read_name: fields[0].to_string(), // Use original read name
                ref_start_pos: alt_hit.pos - 1,
                is_reverse: !alt_hit.strand,
                cigar_str: alt_hit.cigar,
                mapping_quality: 255, // Default MQ for alternative alignments
            };
            process_alignment(alt_alignment, &path_index, &mut stdout)?;
        }
    }

    Ok(())
}

// Two different ways to query a range of steps in a path
fn path_range_cmd(
    path_index: PathIndex,
    path_name: String,
    start: usize,
    end: usize,
) -> Result<()> {
    // Get path info from the indices
    let path = path_index
        .path_names
        .get(&path_name)
        .expect("Path not found");
    let offsets = path_index.path_step_offsets.get(*path).unwrap();
    let steps = path_index.path_steps.get(*path).unwrap();

    // Calculate ranks and cardinality for the range
    let start_rank = offsets.rank(start as u32);
    let end_rank = offsets.rank((end - 1) as u32);
    let cardinality =
        offsets.range_cardinality((start as u32)..((end - 1) as u32));

    println!("start_rank: {start_rank}");
    println!("end_rank: {end_rank}");
    println!("cardinality: {cardinality}");

    println!("------");
    // Calculate how many steps to skip and take
    let skip = (start_rank as usize).checked_sub(1).unwrap_or_default();
    let take = end_rank as usize - skip;

    // Iterate through the steps directly
    println!("step_ix\tnode\tpos");
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

    let pos_range = (start as u32)..((end - 1) as u32);
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
