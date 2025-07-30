# homopolymer_rs
Blazingly fast homopolymer counter for fasta files.

## Requirements
- Linux OS (Ubuntu 24.04.2)
- Rust >= 1.88.0

## Installation
Clone the repository or download the source code. Enter the homopolymer_rs directory and run:<br>
`cargo build --release`

The generated binary is available in `target/release/homopolymer_rs`.

## Usage
Run with:<br>
`homopolymer_rs --fasta <sequences.fasta> --outfile <out.tsv>`

Optional arguments:
<pre>
<b>--min-hp-length</b> [5] - Min homopolymer length to consider.

<b>--strict</b> [false] - Treat uppercase and lowercase nucleotides as different.
</pre>
