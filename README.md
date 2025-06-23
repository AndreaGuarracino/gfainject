# gfainject

`gfainject` is a tool for mapping alignments from `SAM`, `BAM`, `PAF`, or `GBAM` format files to pangenome graphs in a `GFA` (Graphical Fragment Assembly) format graph. It outputs a `GAF` (Graph Alignment Format) file with one record per alignment, mapping each alignment to a sequence of steps in the graph.

## Installation

```shell
git clone https://github.com/AndreaGuarracino/gfainject
cd gfainject
cargo build --release
```

## Usage

```shell
# Convert SAM to GAF
gfainject --gfa ref.gfa --sam aligned.sam > output.gaf

# GFA files can be gzip/bgzip compressed
gfainject --gfa ref.gfa.gz --sam aligned.sam > output.gaf

# SAM input can be piped
samtools view aligned.bam | gfainject --gfa graph.gfa --sam - --alt-hits 5 > output.gaf

# Convert BAM to GAF
gfainject --gfa ref.gfa --bam aligned.bam > output.gaf

# Convert PAF to GAF  
gfainject --gfa ref.gfa --paf aligned.paf > output.gaf

# Convert GBAM to GAF  
gfainject --gfa ref.gfa --gbam aligned.gbam > output.gaf

# Include up to 5 alternative alignments (from XA tag)
gfainject --gfa ref.gfa --bam aligned.bam --alt-hits 5 > output.gaf

# Query a specific range in a path
gfainject --gfa ref.gfa --range "chr1:1000-2000"
```

## Output
`GAF` (Graph Alignment Format) with 12 standard columns plus CIGAR string in `cg:Z` tag.

## Important Notes

- The alignment reference names in the `BAM`, `PAF`, or `GBAM` file must match the path names in the `GFA` file.

## Features

- Fast conversion from `SAM`, `BAM`/`PAF`/`GBAM` to `GAF` format
- Support for compressed `GFA` files (gzip/bgzip)
- Support for alternative alignments via `XA` tags
- Handles paired-end reads with `/1` and `/2` suffixes
- Path range queries for debugging
