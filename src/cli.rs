use clap::{ArgAction, Parser, Subcommand};

#[derive(Parser)]
#[command(
    author = "Sam Dougan",
    version = "0.0.1",
    about = "A kmer-based probe design tool (primary usage: provide fasta, gff, and genes directly)",
    long_about = "This tool is primarily used by providing input files and gene lists directly. \
    The main functionality works with -f/--fasta, -a/--annotation, and -g/--genes arguments. \
    Additional commands like 'build' are available for specific operations.",
    after_help = "\
EXAMPLES:
  GFF mode:
    visiogen gff -f input.fa -a annotation.gff -g gene1,gene2,gene3 -k 50
  
  Other commands:
    visiogen build -i fasta_dir
    visiogen graph -g graph.gfa -t 0.95",
    subcommand_required = true,
    arg_required_else_help = true
)]

pub struct Args {
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 0,
        global = true,
        help = "Number of threads to use for operations (0 = all available cores)"
    )]
    pub threads: usize,

    #[arg(
        short = 'i',
        long = "off_target_directory",
        global = true,
        help = "Directory containing off-target FASTA/index files"
    )]
    pub off_target_directory: Option<String>,

    #[arg(
        long = "max_hits",
        default_value_t = 5,
        global = true,
        help = "Maximum number of index hits to report per kmer"
    )]
    pub max_hits: usize,

    #[arg(
        short = 'r',
        long = "recursive",
        default_value_t = false,
        global = true,
        help = "Recursively search directories for files"
    )]
    pub recursive: bool,

    #[command(flatten)]
    pub kmer_options: KmerOptions,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    Gff(GffArgs),

    Graph(GraphArgs),

    Build(BuildArgs),
}

#[derive(Parser, Debug, Clone)]
pub struct GffArgs {
    #[arg(short = 'a', long = "annotation")]
    pub in_gff: String,

    #[arg(short = 'f', long = "fasta")]
    pub in_fasta: String,

    /// Comma-separated list of gene names
    #[arg(
        short = 'g',
        long = "genes",
        value_delimiter = ',',
        required = true,
        help = "List of gene identifiers comma seperated"
    )]
    pub genes: Vec<String>,
}

#[derive(Parser, Debug, Clone)]
pub struct GraphArgs {
    #[arg(short = 'g', long = "gfa", help = "graph to generate probes from")]
    pub gfa_path: String,
}

#[derive(Parser, Clone)]
pub struct KmerOptions {
    #[arg(
        short = 'k',
        long = "kmer_size",
        default_value_t = 50,
        help = "size of kmer"
    )]
    pub kmer_size: usize,

    #[arg(
        short = 'b',
        long = "center_base",
        help = "Center base leave blank to not consider a center_base"
    )]
    pub center_base: Option<char>,

    #[arg(
        short = 'l',
        long = "min_gc",
        default_value_t = 44,
        help = "Minimum GC content"
    )]
    pub min_gc: usize,

    #[arg(
        short = 'm',
        long = "max_gc",
        default_value_t = 72,
        help = "Maximum GC content"
    )]
    pub max_gc: usize,

    #[arg(
        long = "allow_outside",
        action = ArgAction::SetFalse,
        help = "allow kmers that appear outside target gene"
    )]
    pub allow_outside: bool,

    #[arg(
        long = "skip_gc",
        action = ArgAction::SetTrue,
        help = "skip GC filtering"
    )]
    pub skip_gc: bool,
}

#[derive(Parser, Debug, Clone)]
pub struct BuildArgs {
    /// Use canonical kmers (on by default)
    #[arg(
        short = 'c',
        long = "canonical",
        default_value_t = true,
        help = "Use canonical kmers (default: true)"
    )]
    pub canonical: bool,
}

pub fn parse_args() -> Args {
    let args = Args::try_parse();

    match args {
        Ok(args) => args,
        Err(e) => {
            e.exit();
        }
    }
}
