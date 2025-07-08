use log::*;
use std::collections::HashMap;

use crate::seq;
use crate::GeneKmers;
use crate::Probes;

pub fn generate_gene_kmers(
    genes: &Vec<String>,
    unfiltered_kmers: Vec<Probes>,
    in_gff: &String,
    allow_outside: bool,
) -> Vec<GeneKmers> {
    let results: Vec<GeneKmers> = genes
        .iter()
        .filter_map(|gene| {
            if let Some((start, end, strand)) = seq::coords_from_gene_name(in_gff.clone(), gene) {
                let kmers =
                    seq::filter_hashmap(unfiltered_kmers.clone(), start, end, allow_outside);

                if kmers.is_empty() {
                    info!(
                        "Error: No k-mers uniquely found in gene '{}' {}:{}",
                        gene, start, end
                    );
                    None
                } else {
                    Some(GeneKmers {
                        gene: gene.clone(),
                        start,
                        end,
                        kmers: kmers,
                        strand: strand.to_string(),
                        kmer_hits: HashMap::new(),
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
