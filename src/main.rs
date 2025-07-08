#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use chrono::Local;
use log::*;
use simplelog::*;
use std::fs::File;
use std::path::Path;

use crate::cli::GffArgs;
use crate::cli::KmerOptions;
use crate::index::query_kmers_across_indexes;
use crate::probes::{GeneKmers, Probes};

pub mod cli;
pub mod graph;
pub mod index;
pub mod kmer;
pub mod probes;
pub mod seq;
pub mod utils;

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

fn run_gene_mode(kmer_options: &KmerOptions, gff_args: &GffArgs) -> Vec<GeneKmers> {
    let region = utils::parse_fasta(gff_args.in_fasta.clone()).expect("Failed to parse FASTA file");

    let unfiltered_kmers = Probes::generate_probes(&region, kmer_options.kmer_size, 0); //start at 0 as generating all probes

    info!(
        "total raw kmers: {}",
        unfiltered_kmers.len() / kmer_options.kmer_size
    );

    let gene_kmers = kmer::generate_gene_kmers(
        &gff_args.genes,
        unfiltered_kmers,
        &gff_args.in_gff,
        kmer_options.allow_outside,
    );

    let filtered_all: Vec<GeneKmers> = gene_kmers
        .iter()
        .map(|gk| {
            gk.filter_kmers(
                kmer_options.center_base,
                kmer_options.min_gc,
                kmer_options.max_gc,
                kmer_options.skip_gc,
            )
        })
        .collect();

    filtered_all
}

fn handle_kmer_output(all_kmers: Vec<GeneKmers>, args: &cli::Args, filename_prefix: &str) {
    let kmers_to_write = match args.off_target_directory {
        Some(ref off_target_dir) => {
            match query_kmers_across_indexes(
                Path::new(off_target_dir),
                all_kmers.clone(),
                args.threads,
                args.max_hits,
                args.recursive,
            ) {
                Ok(kmers) => kmers,
                Err(e) => {
                    error!("Failed to query off-target indexes: {}", e);
                    std::process::exit(1);
                }
            }
        }
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
}

fn main() {
    let args = cli::parse_args();

    set_up_logging();

    match &args.command {
        cli::Commands::Gff(gff_args) => {
            let filtered_all = run_gene_mode(&args.kmer_options, gff_args);
            handle_kmer_output(filtered_all, &args, "gff_probes");
        }
        cli::Commands::Graph(graph_args) => {
            let segment_kmers = graph::run_graph_mode(graph_args, args.kmer_options.kmer_size);

            let filtered_all: Vec<GeneKmers> = segment_kmers
                .iter()
                .map(|gk| {
                    gk.filter_kmers(
                        args.kmer_options.center_base,
                        args.kmer_options.min_gc,
                        args.kmer_options.max_gc,
                        args.kmer_options.skip_gc,
                    )
                })
                .collect();

            handle_kmer_output(filtered_all, &args, "graph_probes");
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
    use visiogen::*;

    use super::*;

    #[test]
    fn test_generate_probes() {
        let seq = "ACGTACGT";
        let kmer_size = 3;
        let mut expected: ProbeSet;
        expected.insert("ACG".to_string(), vec![0, 4]);
        expected.insert("CGT".to_string(), vec![1, 5]);
        expected.insert("GTA".to_string(), vec![2]);
        expected.insert("TAC".to_string(), vec![3]);

        let result = Probes::generate_probes(seq, kmer_size, 0);

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
