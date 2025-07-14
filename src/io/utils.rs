use log::*;
use rayon::ThreadPoolBuilder;
use std::path::Path;
use std::path::PathBuf;
use std::sync::Once;
use walkdir::WalkDir;

/// Find all files with the given extensions (e.g., ["fa", "fasta"]) in a directory.
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
