use bio::io::fasta;
use log::*;
use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Once;
use walkdir::WalkDir;

pub fn parse_fasta(fasta_path: String) -> Result<String, io::Error> {
    let read_fasta = isfile(fasta_path)?;
    let fasta_reader = fasta::Reader::new(BufReader::new(read_fasta));

    let mut seq_count = 0;
    let mut result_seq = String::new();

    for record in fasta_reader.records() {
        let record = record?; // incase it errors
        let seq = record.seq();
        let seq_str = String::from_utf8_lossy(seq); // Convert &[u8] to &str

        if seq_count == 0 {
            result_seq = seq_str.to_string();
        }

        seq_count += 1;

        if seq_count > 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Multiple sequences found - currently unsupported",
            ));
        }
    }

    if seq_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Sequence not found",
        ));
    }

    Ok(result_seq)
}

/// Opens a file and returns a more descriptive error if it fails
pub fn isfile(file_path: String) -> io::Result<File> {
    let path = Path::new(&file_path);
    File::open(path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to open file '{}': {}", file_path, e),
        )
    })
}

/// Find all files with the given extensions (e.g., ["fa", "fasta"]) in a directory.
/// If `recursive` is true, searches subdirectories as well.
pub fn find_files_with_extensions(
    directory: &Path,
    extensions: &[&str],
    recursive: bool,
) -> Result<Vec<PathBuf>, Box<dyn std::error::Error>> {
    let ext_set: std::collections::HashSet<_> =
        extensions.iter().map(|e| e.to_lowercase()).collect();

    let files = if recursive {
        WalkDir::new(directory)
            .into_iter()
            .filter_map(Result::ok)
            .filter(|e| e.path().is_file())
            .filter(|e| {
                e.path()
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|ext| ext_set.contains(&ext.to_lowercase()))
                    .unwrap_or(false)
            })
            .map(|e| e.into_path())
            .collect()
    } else {
        std::fs::read_dir(directory)?
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|p| {
                p.is_file()
                    && p.extension()
                        .and_then(|s| s.to_str())
                        .map(|ext| ext_set.contains(&ext.to_lowercase()))
                        .unwrap_or(false)
            })
            .collect()
    };

    Ok(files)
}

pub fn configure_thread_pool(build_threads: usize) {
    static INIT: Once = Once::new();

    INIT.call_once(|| {
        if build_threads > 0 {
            ThreadPoolBuilder::new()
                .num_threads(build_threads)
                .build_global()
                .expect("Failed to build global Rayon thread pool");
        } else {
            info!("Using default number of threads (all logical CPUs)");
            num_cpus::get();
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let sequence = "ATCG";
        let result = reverse_complement(sequence);
        assert_eq!(result, "CGAT");

        let sequence = "AAGCTT";
        let result = reverse_complement(sequence);
        assert_eq!(result, "AAGCTT");

        let sequence = "GGCC";
        let result = reverse_complement(sequence);
        assert_eq!(result, "GGCC");
    }

    #[test]
    fn test_calculate_gc() {
        let sequence = "ATGC";
        let result = calculate_gc(sequence);
        assert_eq!(result, 50); // GC content is 50%

        let sequence = "GGGG";
        let result = calculate_gc(sequence);
        assert_eq!(result, 100); // GC content is 100%

        let sequence = "ATAT";
        let result = calculate_gc(sequence);
        assert_eq!(result, 0); // GC content is 0%
    }

    #[test]
    fn test_gc_content_on_each_half() {
        let kmer = "ATGTCAT";
        let result = gc_content_on_each_half(kmer, kmer.len());
        assert_eq!(result, (33, 33)); // First half: ATG -> 1 G = 33%, Second half: CAT -> 1 G/C = 33%

        let kmer = "GGCCGG";
        let result = gc_content_on_each_half(kmer, kmer.len());
        assert_eq!(result, (100, 100)); // First half and second half both have only GC bases

        let kmer = "ATATAT";
        let result = gc_content_on_each_half(kmer, kmer.len());
        assert_eq!(result, (0, 0)); // Both halves have 0% GC content
    }
}
