use log;
use log::info;
use log::warn;
use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;
/// Ensure we have the proper file ending.
fn check_fasta(fasta: &PathBuf) -> &Path {
    let f_path = Path::new(fasta);

    match f_path.extension().unwrap().to_str().unwrap() {
        "fasta" | "fna" | "fa" => return f_path,
        _ => panic!("File must have a valid fasta extension (cannot be gzipped)."),
    }
}

fn u8_to_char(nt: &u8) -> Option<char> {
    match nt {
        // Default nucleotides.
        b'A' => return Some('A'),
        b'T' => return Some('T'),
        b'C' => return Some('C'),
        b'G' => return Some('G'),
        // Soft masked nucleotides.
        b'a' => return Some('a'),
        b't' => return Some('t'),
        b'c' => return Some('c'),
        b'g' => return Some('g'),
        _ => return None,
    }
}

#[inline]
fn valid_homopolymer(i: usize, j: usize, nt: &u8, min_hp_len: usize, strict: bool) -> bool {
    let valid_len = j - i >= min_hp_len;

    match (strict, u8_to_char(&nt)) {
        (false, _) => return valid_len,
        (true, Some(_)) => return valid_len,
        (true, None) => return false,
    }
}

#[inline]
pub fn find_homopolymers_in_record(
    record: &SequenceRecord,
    min_hp_len: usize,
    strict: bool,
    bufwriter: &mut BufWriter<&File>,
) {
    // Extract sequence information.
    let seq_name = record.id();
    let seq = record.seq();
    let seq_len = seq.len();

    // Skip sequence if shorter than hp len.
    if seq_len < min_hp_len {
        warn!(
            "Skipping record {} (too short)",
            std::str::from_utf8(seq_name).unwrap()
        );
        return;
    }

    let mut i = 0;
    let mut j = 1;

    while i <= seq_len - min_hp_len {
        while j < seq_len && seq[j] == seq[i] {
            j += 1;
        }

        // We have a homopolymer of required length.
        if valid_homopolymer(i, j, &seq[i], min_hp_len, strict) {
            let s = format!(
                "{}\t{}\t{}\t{}\t{}\n",
                std::str::from_utf8(seq_name).unwrap(),
                i,
                j,
                j - i,
                char::from(seq[i])
            );

            bufwriter.write_all(s.as_bytes()).unwrap();
        }

        i = j;
        j += 1;
    }
}

pub fn find_homopolymers_in_fasta(
    fasta: &PathBuf,
    min_hp_len: usize,
    strict: bool,
    outfile: &PathBuf,
) {
    let f_path = check_fasta(fasta);

    // fast(a/q) reader.
    let mut reader = parse_fastx_file(&f_path).unwrap();

    // Output file writer.
    let writer = File::create(outfile).expect("Failed to create output file.");
    let mut bufwriter = BufWriter::new(&writer);

    // Write tsv header.
    let s = format!(
        "{}\t{}\t{}\t{}\t{}\n",
        "contig", "start", "end", "len", "nt"
    );
    bufwriter.write_all(s.as_bytes()).unwrap();

    info!("Finding homopolymers...");
    while let Some(record) = reader.next() {
        match record {
            Ok(record) => find_homopolymers_in_record(&record, min_hp_len, strict, &mut bufwriter),
            Err(e) => {
                warn!("Failed to parse record: {:?}", e);
            }
        }
    }
}
