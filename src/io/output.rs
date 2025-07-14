use crate::cli::Args;
use crate::core::probes::GeneKmers;
use crate::error::{Result, VisiogenError};
use crate::processing::index::query_kmers_across_indexes;
use chrono::Local;
use log::info;
use std::path::Path;

pub fn write_filtered_kmers(
    all_kmers: Vec<GeneKmers>,
    args: &Args,
    filename_prefix: &str,
) -> Result<()> {
    let kmers_to_write = match &args.off_target_directory {
        Some(off_target_dir) => query_kmers_across_indexes(
            Path::new(off_target_dir),
            all_kmers.clone(),
            args.threads,
            args.max_hits,
            args.recursive,
        )
        .map_err(|e| VisiogenError::IndexQueryError(e.to_string()))?,
        None => {
            info!("Skipping off-target check as no off-target directory was provided.");
            all_kmers.clone()
        }
    };

    let timestamp = Local::now().format("%d-%m-%H-%M").to_string();
    let filename = format!("{}_{}.fasta", filename_prefix, timestamp);

    kmers_to_write
        .iter()
        .for_each(|gk| gk.log_and_write_kmers(args.kmer_options.kmer_size, filename.clone()));

    Ok(())
}
