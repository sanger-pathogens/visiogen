use log::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use visiogen::FilteredKmers;

use crate::cli::GraphArgs;
use crate::kmer;

pub fn run_graph_mode(graph_args: &GraphArgs, kmer_size: usize) -> Vec<FilteredKmers> {
    let file = File::open(&graph_args.gfa_path).expect("Failed to open input GFA");
    let reader = BufReader::new(file);

    let mut max_sr = 0usize;
    let mut segment_support: HashMap<String, usize> = HashMap::new();
    let mut segment_seqs: HashMap<String, String> = HashMap::new();
    let mut segment_starts: HashMap<String, u64> = HashMap::new();

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            let seg_id = parts[1].to_string();
            let seq = parts[2].to_string();

            let mut sr_count = 0;
            let mut so_start = 0u64;

            for field in &parts[3..] {
                if field.starts_with("SR:i:") {
                    sr_count = field[5..].parse().unwrap_or(0);
                } else if field.starts_with("SO:i:") {
                    so_start = field[5..].parse().unwrap_or(0);
                }
            }

            max_sr = max_sr.max(sr_count);
            segment_support.insert(seg_id.clone(), sr_count);
            segment_seqs.insert(seg_id.clone(), seq);
            segment_starts.insert(seg_id.clone(), so_start);
        }
    }

    let strain_threshold = ((max_sr as f64) * graph_args.threshold).ceil() as usize;
    info!(
        "Found {max_sr} strains; retaining segments with SR:i ≥ {strain_threshold} (threshold = {:.2})",
        graph_args.threshold
    );

    info!("Total segments available: {}", segment_support.len());

    let core_segments: HashSet<String> = segment_support
        .iter()
        .filter_map(|(seg, &count)| {
            if count >= strain_threshold {
                Some(seg.clone())
            } else {
                None
            }
        })
        .collect();

    info!("Segments passing strain threshold: {}", core_segments.len());

    // Tile each segment individually and produce FilteredKmers per segment
    let filtered_kmers: Vec<FilteredKmers> = core_segments
        .par_iter()
        .filter_map(|seg_id| {
            let seq = segment_seqs.get(seg_id)?;
            let start = *segment_starts.get(seg_id)?;
            let end = start + seq.len() as u64;
            let strand = "+".to_string(); // Not sure just setting strand to +

            let tiled = kmer::tile_segment(seq, start as usize, kmer_size);

            Some(FilteredKmers {
                gene: seg_id.clone(),
                start,
                end,
                kmers: tiled,
                strand,
                kmer_hits: HashMap::new(),
            })
        })
        .collect();

    let total_kmers: usize = filtered_kmers.iter().map(|f| f.kmers.len()).sum();
    info!(
        "Generated kmers for {} segments (total raw kmers: {}, avg per segment: {:.2})",
        filtered_kmers.len(),
        total_kmers,
        total_kmers as f64 / filtered_kmers.len().max(1) as f64
    );

    filtered_kmers
}
