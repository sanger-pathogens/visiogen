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

    let mut segment_seqs: HashMap<String, String> = HashMap::new();
    let mut segment_starts: HashMap<String, u64> = HashMap::new();
    let mut segment_connectivity: HashMap<String, SegmentConnectivity> = HashMap::new();

    // Parse segments and links
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        let parts: Vec<&str> = line.split('\t').collect();

        if line.starts_with('S') {
            let seg_id = parts[1].to_string();
            let seq = parts[2].to_string();
            let mut so = 0u64;

            for field in &parts[3..] {
                if field.starts_with("SO:i:") {
                    so = field[5..].parse().unwrap_or(0);
                }
            }

            segment_seqs.insert(seg_id.clone(), seq);
            segment_starts.insert(seg_id.clone(), so);
        } else if line.starts_with('L') {
            let from = parts[1].to_string();
            let to = parts[3].to_string();

            segment_connectivity.entry(from.clone())
                .or_default()
                .outgoing
                .insert(to.clone());
            segment_connectivity.entry(to.clone())
                .or_default()
                .incoming
                .insert(from.clone());
        }
    }

    // Start walking from "s1"
    let mut current = "s1".to_string();
    let mut visited = HashSet::new();
    let mut unique_path_segments = Vec::new();

    loop {
        if visited.contains(&current) {
            break;
        }
        visited.insert(current.clone());

        let conn = segment_connectivity.get(&current);
        if conn.is_none() {
            break;
        }

        let conn = conn.unwrap();
        let outgoing = &conn.outgoing;
        let incoming = &conn.incoming;

        // Unique if exactly one in and one out
        if incoming.len() <= 1 && outgoing.len() == 1 {
            unique_path_segments.push(current.clone());
        }

        // Choose the lowest-numbered next segment if multiple
        let next = outgoing.iter()
            .min_by_key(|s| s[1..].parse::<u64>().unwrap_or(u64::MAX)); // strip 's' and parse

        match next {
            Some(n) => current = n.clone(),
            None => break,
        }
    }

    info!("Identified {} unique path segments", unique_path_segments.len());

    // Generate kmers for each segment on the path
    let filtered_kmers: Vec<FilteredKmers> = unique_path_segments
        .par_iter()
        .filter_map(|seg_id| {
            let seq = segment_seqs.get(seg_id)?;
            let start = *segment_starts.get(seg_id)?;
            let end = start + seq.len() as u64;

            Some(FilteredKmers {
                gene: seg_id.clone(),
                start,
                end,
                kmers: kmer::tile_segment(seq, start as usize, kmer_size),
                strand: "+".to_string(),
                kmer_hits: HashMap::new(),
            })
        })
        .collect();

    let total_kmers: usize = filtered_kmers.iter().map(|f| f.kmers.len()).sum();
    info!(
        "Generated kmers for {} segments (total kmers: {}, avg per segment: {:.2})",
        filtered_kmers.len(),
        total_kmers,
        total_kmers as f64 / filtered_kmers.len().max(1) as f64
    );

    filtered_kmers
}