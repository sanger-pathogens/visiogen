#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

mod cli;
mod core;
mod error;
mod io;
mod logging;
mod processing;

use crate::cli::{parse_args, Args, BuildArgs, Commands, GffArgs, KmerOptions};
use crate::core::probes::{GeneKmers, Probes};
use crate::error::{Result, VisiogenError};
use crate::io::output;
use crate::processing::{gff, graph, index};
use log::info;
use std::collections::HashMap;

fn run(args: Args) -> Result<()> {
    match &args.command {
        Commands::Gff(gff_args) => run_probe_command(&args, gff_args),
        Commands::Build(build_args) => run_build_command(&args, build_args),
    }
}

fn run_probe_command(args: &Args, gff_args: &GffArgs) -> Result<()> {
    let _gene_coords =
        gff::coords_from_gene_name(&gff_args.in_gff, &gff_args.genes).map_err(|e| {
            VisiogenError::GffParseError(format!(
                "off_target_directory required for build command {}",
                e
            ))
        })?;

    let graph = graph::parse_gfa_file(&args.gfa_path)
        .map_err(|e| VisiogenError::GfaParseError(format!("Failed to read GFA file: {}", e)))?;

    let segment_kmers: Vec<GeneKmers> = graph
        .core_segment_structs()
        .iter()
        .map(|segment| GeneKmers {
            gene: segment.name.clone(),
            start: 1,
            end: 1 + segment.sequence.len() as u64,
            kmers: Probes::generate_probes(&segment.sequence, args.kmer_options.kmer_size, 0),
            strand: "+".to_string(),
            kmer_hits: HashMap::new(),
        })
        .collect();

    let total_kmers: usize = segment_kmers.iter().map(|f| f.kmers.len()).sum();
    info!(
        "Generated kmers for {} segments (total kmers: {}, avg per segment: {:.2})",
        segment_kmers.len(),
        total_kmers,
        total_kmers as f64 / segment_kmers.len().max(1) as f64
    );

    let filtered_kmers = apply_kmer_filters(segment_kmers, &args.kmer_options);

    let final_probes = select_best_probes(filtered_kmers, args.n_count);

    output::write_filtered_kmers(final_probes, args, "probes")?;

    Ok(())
}

fn run_build_command(args: &Args, build_args: &BuildArgs) -> Result<()> {
    let off_target_dir = args.off_target_directory.as_ref().ok_or_else(|| {
        VisiogenError::MissingArgument(
            "off_target_directory required for build command".to_string(),
        )
    })?;

    index::build_indexes_for_all_fastas(
        std::path::Path::new(off_target_dir),
        args.threads,
        build_args.canonical,
        args.recursive,
    )
    .map_err(|e| {
        VisiogenError::IndexBuildError(format!("Failed to build indexes for fastas {}", e))
    })?;

    Ok(())
}

fn apply_kmer_filters(gene_kmers: Vec<GeneKmers>, kmer_options: &KmerOptions) -> Vec<GeneKmers> {
    gene_kmers
        .iter()
        .map(|gk| {
            gk.filter_kmers(
                kmer_options.center_base,
                kmer_options.min_gc,
                kmer_options.max_gc,
                kmer_options.skip_gc,
            )
        })
        .collect()
}

fn select_best_probes(gene_kmers: Vec<GeneKmers>, n_count: u16) -> Vec<GeneKmers> {
    gene_kmers
        .iter()
        .map(|gk| gk.best_probes(n_count))
        .collect()
}

fn main() {
    let args = parse_args();
    logging::set_up_logging();

    if let Err(e) = run(args) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
