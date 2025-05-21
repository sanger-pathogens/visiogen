use bio_types::strand::Strand;
use dashmap::DashMap;
use kmer_visium::FilteredKmers;
use log::*;
use rayon::prelude::*;
use std::collections::HashMap;

use crate::seq;

pub fn tile_string(seq: &str, kmer_size: usize) -> HashMap<String, Vec<usize>> {
    let kmers = DashMap::new();
    let seq_len = seq.len();

    //pool to parralise
    let indices: Vec<usize> = (0..=seq_len - kmer_size).collect();

    indices.into_par_iter().for_each(|i| {
        let kmer = &seq[i..i + kmer_size];
        kmers
            .entry(kmer.to_string())
            .or_insert_with(Vec::new)
            .push(i);
    });

    //swap the dashmap back into the final HashMap
    kmers.into_iter().collect()
}

pub fn tile_segment(
    seq: &str,
    start_offset: usize,
    kmer_size: usize,
) -> HashMap<String, Vec<usize>> {
    let kmers = DashMap::new();
    let seq_len = seq.len();

    if seq_len < kmer_size {
        return HashMap::new(); // skip too-short segments
    }

    let indices: Vec<usize> = (0..=seq_len - kmer_size).collect();

    indices.into_par_iter().for_each(|i| {
        let kmer = &seq[i..i + kmer_size];
        kmers
            .entry(kmer.to_string())
            .or_insert_with(Vec::new)
            .push(i + start_offset);
    });

    kmers.into_iter().collect()
}

pub fn generate_gene_kmers(
    genes: Vec<String>,
    unfiltered_kmers: HashMap<String, Vec<usize>>,
    in_gff: String,
    allow_outside: bool,
) -> Vec<FilteredKmers> {
    let results: Vec<FilteredKmers> = genes
        .iter()
        .filter_map(|gene| {
            if let Some((start, end, strand)) = seq::coords_from_gene_name(in_gff.clone(), gene) {
                let kmers = seq::filter_hashmap(&unfiltered_kmers, start, end, allow_outside);

                if kmers.is_empty() {
                    info!(
                        "Error: No k-mers uniquely found in gene '{}' {}:{}",
                        gene, start, end
                    );
                    None
                } else {
                    // If the strand is negative, apply reverse_complement to the kmers
                    let processed_kmers = match strand {
                        Strand::Reverse => kmers
                            .into_iter()
                            .map(|(kmer, count)| (seq::reverse_complement(&kmer), count))
                            .collect(),
                        _ => kmers,
                    };
                    Some(FilteredKmers {
                        gene: gene.clone(),
                        start,
                        end,
                        kmers: processed_kmers,
                        strand: strand.to_string(),
                    })
                }
            } else {
                info!("Gene {} not found", gene);
                None
            }
        })
        .collect();
    results
}
