#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use chrono::Local;
use log::*;
use rayon::prelude::*;
use simplelog::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::cli::GffArgs;
use crate::cli::KmerOptions;
use crate::index::query_kmers_across_indexes;
use visiogen::FilteredKmers;

mod cli;
mod graph;
mod index;
mod kmer;
mod seq;
mod utils;

fn set_up_logging() {
    let current_time = Local::now().format("%m-%d_%H-%M-%S").to_string();
    let log_filename = format!("visiogen_{}.log", current_time);

    CombinedLogger::init(vec![
        TermLogger::new(
            LevelFilter::Warn,
            Config::default(),
            TerminalMode::Mixed,
            ColorChoice::Auto,
        ),
        WriteLogger::new(
            LevelFilter::Info,
            Config::default(),
            File::create(log_filename).unwrap(),
        ),
    ])
    .unwrap();
}

fn filter_kmers(
    filtered_kmers_vec: &Vec<visiogen::FilteredKmers>,
    kmer_size: usize,
    center_base: Option<char>,
    min_gc: usize,
    max_gc: usize,
    skip_gc: bool,
) -> Vec<visiogen::FilteredKmers> {
    filtered_kmers_vec
        .iter()
        .map(|filtered_kmers| {
            let valid_kmers: HashMap<String, Vec<usize>> = filtered_kmers
                .kmers
                .par_iter()
                .filter_map(|(kmer, positions)| {
                    if let Some(center) = center_base {
                        if kmer.chars().nth(kmer_size / 2 - 1) != Some(center) {
                            debug! {"kmer: {} failed as middle base was not {}", kmer, center};
                            return None;
                        }
                    }

                    if !skip_gc {
                        let (gc_first, gc_second) = seq::gc_content_on_each_half(kmer, kmer_size);

                        if !(gc_first >= min_gc
                            && gc_first <= max_gc
                            && gc_second >= min_gc
                            && gc_second <= max_gc)
                        {
                            debug! {"kmer: {} failed as gc {} - {} - {} was out of range {} - {}",
                            kmer, gc_first, center_base.unwrap_or('-'), gc_second, min_gc, max_gc}; // Use '-' as placeholder if no center_base
                            return None;
                        }
                    }

                    Some((kmer.clone(), positions.clone()))
                })
                .collect();

            visiogen::FilteredKmers {
                gene: filtered_kmers.gene.clone(),
                start: filtered_kmers.start.clone(),
                end: filtered_kmers.end.clone(),
                kmers: valid_kmers,
                strand: filtered_kmers.strand.clone(),
                kmer_hits: HashMap::new(),
            }
        })
        .collect()
}

fn log_kmers_with_coords(filtered_kmers: &visiogen::FilteredKmers, kmer_size: usize) {
    for (kmer, start_coords) in &filtered_kmers.kmers {
        for &start in start_coords {
            let end = if filtered_kmers.strand == "-" {
                // If strand is '-', go backwards
                if start >= kmer_size {
                    start - kmer_size
                } else {
                    0 // Ensure end coordinate doesn't go negative
                }
            } else {
                // Otherwise go forwards
                start + kmer_size
            };

            info!("{},{},{}", kmer, start, end);
        }
    }
}

fn write_all_keys_to_file(all_filtered_kmers: &[FilteredKmers]) {
    let timestamp = Local::now().format("%m-%d_%H-%M-%S").to_string();
    let filename = format!("{}.fasta", timestamp);
    let mut final_file = File::create(&filename).expect("Failed to create FASTA file");

    for filtered in all_filtered_kmers {
        for (i, (kmer, coords)) in filtered.kmers.iter().enumerate() {
            let coords_str = coords
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(",");

            writeln!(
                final_file,
                ">{}_{}    {} : {} copies",
                filtered.gene,
                i + 1,
                coords_str,
                coords.len()
            )
            .expect("Failed to write FASTA header");
            writeln!(final_file, "{}", kmer).expect("Failed to write FASTA sequence");
        }
    }

    info!("Wrote all kmers to file: {}", filename);
}

fn log_and_write_kmers(all_filtered_kmers: &[FilteredKmers], kmer_size: usize) {
    // Write all kmers to a timestamped .fasta file
    write_all_keys_to_file(all_filtered_kmers);

    // Log metadata and optionally detailed coordinates
    for filtered in all_filtered_kmers {
        info!("Gene: {}", filtered.gene);
        info!("Strand: {}", filtered.strand);
        info!("Start: {}", filtered.start);
        info!("End: {}", filtered.end);
        info!("Total: {}", filtered.kmers.len());

        if log_enabled!(Level::Debug) {
            log_kmers_with_coords(filtered, kmer_size);
        }
    }
}

fn search_kmers(kmer_options: &KmerOptions, gff_args: GffArgs) -> Vec<FilteredKmers> {
    let region = utils::parse_fasta(gff_args.in_fasta.clone()).expect("Failed to parse FASTA file");

    let unfiltered_kmers = kmer::tile_string(&region, kmer_options.kmer_size);

    info!(
        "total raw kmers: {}",
        unfiltered_kmers.len() / kmer_options.kmer_size
    );

    let gene_kmers = kmer::generate_gene_kmers(
        gff_args.genes,
        unfiltered_kmers,
        gff_args.in_gff,
        kmer_options.allow_outside,
    );

    let all_filtered_kmers = filter_kmers(
        &gene_kmers,
        kmer_options.kmer_size,
        kmer_options.center_base,
        kmer_options.min_gc,
        kmer_options.max_gc,
        kmer_options.skip_gc,
    );

    all_filtered_kmers
}

fn main() {
    let args = cli::parse_args();

    set_up_logging();

    match args.command {
        cli::Commands::Gff(gff_args) => {
            let all_filtered_kmers = search_kmers(&args.kmer_options, gff_args);

            let kmers_to_write = match args.off_target_directory {
                Some(ref off_target_dir) => {
                    match query_kmers_across_indexes(
                        Path::new(off_target_dir),
                        all_filtered_kmers.clone(),
                        args.threads,
                        args.max_hits,
                        args.recursive,
                    ) {
                        Ok(kmers) => kmers,
                        Err(e) => {
                            log::error!("Failed to query off-target indexes: {}", e);
                            std::process::exit(1);
                        }
                    }
                }
                None => {
                    log::info!(
                        "Skipping off-target check as no off-target directory was provided."
                    );
                    all_filtered_kmers.clone()
                }
            };

            log_and_write_kmers(&kmers_to_write, args.kmer_options.kmer_size);
        }
        cli::Commands::Graph(graph_args) => {
            let segment_kmers = graph::run_graph_mode(&graph_args, args.kmer_options.kmer_size);

            let all_filtered_kmers = filter_kmers(
                &segment_kmers,
                args.kmer_options.kmer_size,
                args.kmer_options.center_base,
                args.kmer_options.min_gc,
                args.kmer_options.max_gc,
                args.kmer_options.skip_gc,
            );

            let kmers_to_write = match args.off_target_directory {
                Some(ref off_target_dir) => {
                    match query_kmers_across_indexes(
                        Path::new(off_target_dir),
                        all_filtered_kmers.clone(),
                        args.threads,
                        args.max_hits,
                        args.recursive,
                    ) {
                        Ok(kmers) => kmers,
                        Err(e) => {
                            log::error!("Failed to query off-target indexes: {}", e);
                            std::process::exit(1);
                        }
                    }
                }
                None => {
                    log::info!(
                        "Skipping off-target check as no off-target directory was provided."
                    );
                    all_filtered_kmers.clone()
                }
            };

            log_and_write_kmers(&kmers_to_write, args.kmer_options.kmer_size);
        }
        cli::Commands::Build(build_args) => {
            info!("indexing: {:#?}", args.off_target_directory);
            index::build_indexes_for_all_fastas(
                Path::new(args.off_target_directory.as_ref().unwrap()),
                args.threads,
                build_args.canonical,
                args.recursive,
            )
            .expect("Failed to build index");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tile_string() {
        let seq = "ACGTACGT";
        let kmer_size = 3;
        let mut expected: HashMap<String, Vec<usize>> = HashMap::new();
        expected.insert("ACG".to_string(), vec![0, 4]);
        expected.insert("CGT".to_string(), vec![1, 5]);
        expected.insert("GTA".to_string(), vec![2]);
        expected.insert("TAC".to_string(), vec![3]);

        let result = kmer::tile_string(seq, kmer_size);

        // Check that all keys in expected exist in result and have the same values (ignoring order)
        assert_eq!(result.len(), expected.len());

        for (key, expected_vec) in expected.iter() {
            assert!(
                result.contains_key(key),
                "Key '{}' not found in result",
                key
            );

            let mut result_vec = result.get(key).unwrap().clone();
            let mut expected_vec = expected_vec.clone();

            // Sort both vectors to ensure order does not matter
            result_vec.sort_unstable();
            expected_vec.sort_unstable();

            assert_eq!(result_vec, expected_vec, "Values for key '{}' differ", key);
        }
    }
}
