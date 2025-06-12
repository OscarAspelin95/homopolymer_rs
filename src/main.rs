use bio::io::fasta::{Reader, Record};
use clap::Parser;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(long_about = "Finds homopolymers in fasta files.")]
struct CommandArgs {
    #[arg(short, long)]
    fasta: PathBuf,

    #[arg(long, default_value = "5",  value_parser= clap::value_parser!(u8).range(5..20))]
    pub min_hp_length: u8,
}

#[inline]
pub fn u8_to_char(x: &u8) -> char {
    match x {
        // A/a.
        b'A' | b'a' => return 'A',
        // T/t.
        b'T' | b't' => return 'T',
        // C/c.
        b'C' | b'c' => return 'C',
        // G/g.
        b'G' | b'g' => return 'G',
        // N only.
        b'N' => return 'N',
        _ => panic!("Invalid ascii nucleotide {x}"),
    };
}

#[inline]
pub fn find_homopolymers<'a>(record: &'a Record, min_hp_len: usize) {
    // Store results.

    // Extract sequence information.
    let seq_name = record.id();
    let seq = record.seq();
    let seq_len = seq.len();

    // Skip sequence if shorter than hp len.
    if seq_len < min_hp_len {
        return;
    }

    let mut i = 0;
    let mut j = 1;

    while i <= seq_len {
        while j < seq_len && seq[j] == seq[i] {
            j += 1;
        }

        if j - i >= min_hp_len {
            // contig name  start   end length  nucleotide
            println!(
                "{}\t{}\t{}\t{}\t{}",
                seq_name,
                i,
                j,
                j - i,
                u8_to_char(&seq[i])
            );
        }

        i = j;
        j += 1;
    }
}

fn check_fasta(args: &CommandArgs) -> &Path {
    let f_path = Path::new(&args.fasta);

    match f_path.extension().unwrap().to_str().unwrap() {
        "fasta" | "fna" | "fa" => return f_path,
        _ => panic!("File must have a valid fasta extension (cannot be gzipped)."),
    }
}
fn main() {
    let args: CommandArgs = CommandArgs::parse();

    let f_path = check_fasta(&args);

    // Extract args.
    let fasta_file = File::open(f_path).unwrap();
    let min_hp_len: usize = args.min_hp_length as usize;

    // Read fasta.
    let bufreader = BufReader::new(fasta_file);
    let reader = Reader::from_bufread(bufreader);
    let mut records = reader.records();

    println!("{}\t{}\t{}\t{}\t{}", "contig", "start", "end", "len", "nt");

    // This is actually not ideal because the loop will terminate prematurely
    // without an error if a faulty record is seen.
    // A better alternative is to loop over all records and call .unwrap().
    while let Some(Ok(record)) = records.next() {
        find_homopolymers(&record, min_hp_len);
    }
}
