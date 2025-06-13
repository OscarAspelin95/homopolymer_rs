use clap::Parser;
use std::path::PathBuf;

mod homopolymers;
use homopolymers::print_homopolymers_in_fasta;

#[derive(Parser, Debug)]
#[command(long_about = "Finds homopolymers in fasta files.")]
struct CommandArgs {
    #[arg(short, long)]
    fasta: PathBuf,

    #[arg(long, default_value = "5",  value_parser= clap::value_parser!(u8).range(5..))]
    min_hp_length: u8,

    #[arg(long, action)]
    strict: bool,
}

/// TODO - add .gz support.
/// TODO - (possibly) add rayon multithreading (either per record or per fasta file).
/// TODO - add a "--strict" mode to only allow/print unambiguous nucleotides ATCGatcg.
fn main() {
    let args: CommandArgs = CommandArgs::parse();

    print_homopolymers_in_fasta(&args.fasta, args.min_hp_length as usize, args.strict);
}
