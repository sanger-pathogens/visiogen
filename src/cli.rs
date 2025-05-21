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
    visiogen build -i fasta_dir -o output.idx
    visiogen graph -g graph.gfa -t 0.95",
    subcommand_required = true,
    arg_required_else_help = true
)]

pub struct Args {
    /// Number of threads to use (0 = use all available cores)
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 0,
        global = true,
        help = "Number of threads to use for operations (0 = all available cores)"
    )]
    pub threads: usize,

    /// Directory containing off-target FASTA/index files
    #[arg(
        short = 'i',
        long = "off_target_directory",
        global = true,
        help = "Directory containing off-target FASTA/index files"
    )]
    pub off_target_directory: Option<String>,

    /// Maximum number of index hits to report per kmer (used for querying)
    #[arg(
        long = "max_hits",
        default_value_t = 5,
        global = true,
        help = "Maximum number of index hits to report per kmer"
    )]
    pub max_hits: usize,

    /// Recursively search directories (used for indexing and querying)
    #[arg(
        short = 'r',
        long = "recursive",
        default_value_t = false,
        global = true,
        help = "Recursively search directories for files"
    )]
    pub recursive: bool,

    /// Kmer options
    #[command(flatten)]
    pub kmer_options: KmerOptions,

    /// Subcommand
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Work with GFF annotations and FASTA files
    Gff(GffArgs),

    /// Graph based args
    Graph(GraphArgs),

    /// Build an off target index
    Build(BuildArgs),
}

#[derive(Parser)]
pub struct GffArgs {
    /// Input file - GFF
    #[arg(short = 'a', long = "annotation")]
    pub in_gff: String, // Changed from Option<String> to String since it's required in GFF mode

    /// Input file - FASTA
    #[arg(short = 'f', long = "fasta")]
    pub in_fasta: String, // Changed from Option<String> to String

    /// Comma-separated list of gene names
    #[arg(short = 'g', long = "genes", value_delimiter = ',', required = true)]
    pub genes: Vec<String>,
}

#[derive(Parser, Debug)]
pub struct GraphArgs {
    /// Input GFA file
    #[arg(short = 'g', long = "gfa")]
    pub gfa_path: String,

    /// Core segment inclusion threshold (fraction)
    #[arg(short = 't', long = "threshold", default_value_t = 0.95)]
    pub threshold: f64,
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

#[derive(Parser, Debug)]
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
