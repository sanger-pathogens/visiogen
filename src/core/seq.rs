use bio::alphabets::dna::revcomp;

use crate::core::probes::ProbeSet;

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
