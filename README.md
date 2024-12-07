# gfainject

`gfainject` is a tool for mapping alignments from BAM, GBAM, or PAF format files to pangenome graphs in a GFA (Graphical Fragment Assembly) format graph. It outputs a GAF (Graph Alignment Format) file with one record per alignment, mapping each alignment to a sequence of steps in the graph.

## Features

- Supports BAM, PAF, and GBAM input formats
- Maps alignments to GFA paths
- Outputs results in GAF format
- Efficient processing of large datasets

## Requirements

- Rust programming environment
- Input files:
  - GFA file containing the pangenome graph
  - BAM, PAF, or GBAM file with alignments

## Usage

### Building the tool

```sh
git clone https://github.com/AndreaGuarracino/gfainject
cd gfainject
cargo build --release
```

### Running gfainject

For BAM input:
```sh
./target/release/gfainject --gfa <path_to_gfa> --bam <path_to_bam>
```

For GBAM input:

```sh
./target/release/gfainject --gfa <path_to_gfa> --gbam <path_to_gbam>
```

For PAF input:
```sh
./target/release/gfainject --gfa <path_to_gfa> --paf <path_to_paf>
```

### Output

The tool generates a GAF file to `stdout`. You can redirect this to a file:

```sh
./target/release/gfainject --gfa input.gfa --bam input.bam > output.gaf
```

## Important Notes

- The alignment reference names in the BAM, PAF, or GBAM file must match the path names in the GFA file.

## Additional Features

- Path range querying: You can query a specific range within a path using the `--range` option with format "path_name:start-end".

```sh
./target/release/gfainject --gfa <path_to_gfa> --range "path_name:1000-2000"
```
