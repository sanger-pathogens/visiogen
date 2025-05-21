use bio::alphabets::dna::revcomp;
use bio::io::gff;
use bio_types::strand::Strand;
use std::collections::HashMap;
use std::io::BufReader;

use crate::utils;

pub fn reverse_complement(sequence: &str) -> String {
    // Get reverse complement and convert it back to String
    revcomp(sequence.as_bytes())
        .iter()
        .map(|&b| b as char)
        .collect()
}

fn calculate_gc(sequence: &str) -> usize {
    let total_length = sequence.len();
    let gc_count = sequence
        .chars()
        .filter(|&c| c == 'G' || c == 'g' || c == 'C' || c == 'c')
        .count();

    let gc_content_percentage = (gc_count * 100) / total_length;

    gc_content_percentage
}

pub fn gc_content_on_each_half(kmer: &str, kmer_size: usize) -> (usize, usize) {
    let mid_index = kmer_size / 2;

    let first_gc = calculate_gc(&kmer[..mid_index]);
    let second_gc = calculate_gc(&kmer[mid_index + 1..]);

    (first_gc, second_gc)
}

pub fn filter_hashmap<'a>(
    kmer_hash: &'a HashMap<String, Vec<usize>>,
    start: u64,
    end: u64,
    allow_outside: bool,
) -> HashMap<String, Vec<usize>> {
    kmer_hash
        .iter() // Iterate over references to the entries
        .filter_map(|(key, indices)| {
            if allow_outside {
                // Use `all` to check if all indices are within range
                if indices
                    .iter()
                    .all(|&index| index >= (start as usize) && index <= (end as usize))
                {
                    Some((key.to_string(), indices.to_vec()))
                } else {
                    None
                }
            } else {
                // Use `any` to check if any index is within range
                if indices
                    .iter()
                    .any(|&index| index >= (start as usize) && index <= (end as usize))
                {
                    Some((key.to_string(), indices.to_vec()))
                } else {
                    None
                }
            }
        })
        .collect()
}

pub fn coords_from_gene_name(gff_path: String, gene: &str) -> Option<(u64, u64, Strand)> {
    let gff_file = utils::isfile(gff_path).ok()?;
    let mut gff_reader = gff::Reader::new(BufReader::new(gff_file), gff::GffType::GFF3);

    for record in gff_reader.records() {
        let rec = record.expect("Error reading record.");
        if let Some(attributes) = rec.attributes().get("Name") {
            if attributes == gene {
                return Some((
                    *rec.start(),
                    *rec.end(),
                    rec.strand().unwrap_or(Strand::Forward),
                ));
            }
        }
    }
    None
}
