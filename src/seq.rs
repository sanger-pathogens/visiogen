use bio::alphabets::dna::revcomp;
use bio::io::gff;
use bio_types::strand::Strand;
use std::io::BufReader;

use crate::probes::ProbeSet;
use crate::utils;

pub fn reverse_complement(sequence: &str) -> String {
    // Get reverse complement and convert it back to String
    revcomp(sequence.as_bytes())
        .iter()
        .map(|&b| b as char)
        .collect()
}

pub fn filter_hashmap(probes: ProbeSet, start: u64, end: u64, allow_outside: bool) -> ProbeSet {
    probes
        .into_iter()
        .filter(|probe| {
            if allow_outside {
                // All locations must be within the range
                probe
                    .locations
                    .iter()
                    .all(|&pos| pos >= start as usize && pos <= end as usize)
            } else {
                // At least one location must be within the range
                probe
                    .locations
                    .iter()
                    .any(|&pos| pos >= start as usize && pos <= end as usize)
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
