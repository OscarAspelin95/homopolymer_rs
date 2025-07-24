use clap::Parser;
use std::path::PathBuf;
mod homopolymers;
use homopolymers::find_homopolymers_in_fasta;

#[derive(Parser, Debug)]
#[command(long_about = "Finds homopolymers in the forward strand of fasta entries.")]
struct CommandArgs {
    #[arg(short, long)]
    fasta: PathBuf,

    #[arg(long, default_value = "5",  value_parser= clap::value_parser!(u8).range(5..))]
    min_hp_length: u8,

    #[arg(long, action)]
    strict: bool,
}

/// TODO - add tests.
fn main() {
    let args: CommandArgs = CommandArgs::parse();

    find_homopolymers_in_fasta(&args.fasta, args.min_hp_length as usize, args.strict);
}
