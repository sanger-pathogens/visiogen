use log::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::io::Write;

pub type ProbeSet = Vec<Probes>;

#[derive(Debug, Clone)]
pub struct GeneKmers {
    pub gene: String,
    pub start: u64,
    pub end: u64,
    pub kmers: ProbeSet, // type ProbeSet = Vec<Probes>
    pub strand: String,
    pub kmer_hits: HashMap<String, Vec<String>>,
}

impl GeneKmers {
    pub fn filter_kmers(
        &self,
        center_base: Option<char>,
        min_gc: usize,
        max_gc: usize,
        skip_gc: bool,
    ) -> GeneKmers {
        let valid_kmers: Vec<Probes> = self
            .kmers
            .par_iter()
            .filter(|probe| {
                let junction_matches = match center_base {
                    Some(base) => probe.junction_base == base,
                    None => true,
                };

                let first_gc_valid = min_gc <= probe.first_half_gc && probe.first_half_gc <= max_gc;
                let second_gc_valid =
                    min_gc <= probe.second_half_gc && probe.second_half_gc <= max_gc;
                let gc_valid = skip_gc || (first_gc_valid && second_gc_valid);

                junction_matches && gc_valid
            })
            .cloned()
            .collect();

        GeneKmers {
            gene: self.gene.clone(),
            start: self.start,
            end: self.end,
            kmers: valid_kmers,
            strand: self.strand.clone(),
            kmer_hits: HashMap::new(),
        }
    }

    pub fn log_kmers_with_coords(&self, kmer_size: usize) {
        for probe in &self.kmers {
            for &start in &probe.locations {
                let end = if self.strand == "-" {
                    if start >= kmer_size {
                        start - kmer_size
                    } else {
                        0
                    }
                } else {
                    start + kmer_size
                };

                info!("{},{},{}", probe.kmer, start, end);
            }
        }
    }

    pub fn write_all_keys_to_file(&self, filename: String) {
        let mut final_file = OpenOptions::new()
            .create(true) // Creates the file if it doesn't exist
            .append(true) // Appends to the file if it does exist
            .open(&filename)
            .expect("Failed to open or create FASTA file");

        for (i, probe) in self.kmers.iter().enumerate() {
            let coords_str = probe
                .locations
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(",");

            writeln!(
                final_file,
                ">{}_{}    {} : {} copies",
                self.gene,
                i + 1,
                coords_str,
                probe.locations.len()
            )
            .expect("Failed to write FASTA header");

            writeln!(final_file, "{}", probe.kmer).expect("Failed to write FASTA sequence");
        }
    }

    pub fn log_and_write_kmers(&self, kmer_size: usize, filename: String) {
        Self::write_all_keys_to_file(self, filename);

        info!(
            "Gene: {}, Strand: {}, Start: {}, End: {}, Total: {}",
            self.gene,
            self.strand,
            self.start,
            self.end,
            self.kmers.len()
        );

        if log_enabled!(Level::Debug) {
            self.log_kmers_with_coords(kmer_size);
        }
    }
}

#[derive(Debug, Clone)]
pub struct Probes {
    pub kmer: String,
    pub locations: Vec<usize>,
    pub first_half_gc: usize,
    pub second_half_gc: usize,
    pub complexity: f64,
    pub junction_base: char,
    pub score: Option<f64>,
}

impl Probes {
    fn new(kmer: String, locations: Vec<usize>) -> Self {
        let first_half_gc = Self::calculate_gc(&kmer[..kmer.len() / 2]);
        let second_half_gc = Self::calculate_gc(&kmer[kmer.len() / 2..]);
        let complexity = Self::score_homopolymer_repeats(&kmer);
        let junction_base = kmer.chars().nth(24).unwrap_or('N');
        let score = None;

        Self {
            kmer,
            locations,
            first_half_gc,
            second_half_gc,
            complexity,
            junction_base,
            score,
        }
    }

    pub fn generate_probes(seq: &str, kmer_size: usize, start_offset: usize) -> ProbeSet {
        let mut kmers: HashMap<String, Vec<usize>> = HashMap::new();

        for i in 0..=seq.len() - kmer_size {
            let kmer = &seq[i..i + kmer_size];
            kmers
                .entry(kmer.to_string())
                .or_default()
                .push(i + start_offset);
        }

        kmers
            .into_iter()
            .map(|(kmer, locations)| Self::new(kmer, locations))
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

    /// Returns a complexity score between 0.0 (very repetitive) and 1.0 (diverse)
    fn score_homopolymer_repeats(seq: &str) -> f64 {
        let mut max_run = 1;
        let mut current_run = 1;
        let mut prev_char = None;

        for c in seq.chars() {
            if Some(c) == prev_char {
                current_run += 1;
                max_run = max_run.max(current_run);
            } else {
                current_run = 1;
            }
            prev_char = Some(c);
        }

        let length = seq.len() as f64;
        let repeat_fraction = max_run as f64 / length;

        1.0 - repeat_fraction
    }
}
