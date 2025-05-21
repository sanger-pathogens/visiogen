use log::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use visiogen::FilteredKmers;

use crate::cli::GraphArgs;
use crate::kmer;

#[derive(Default)]
struct SegmentConnectivity {
    incoming: HashSet<String>,
    outgoing: HashSet<String>,
}

pub fn run_graph_mode(graph_args: &GraphArgs, kmer_size: usize) -> Vec<FilteredKmers> {

    let file = File::open(&graph_args.gfa_path).expect("Failed to open input GFA");
    let reader = BufReader::new(file);

    let mut segment_seqs: HashMap<String, String> = HashMap::new();
    let mut segment_starts: HashMap<String, u64> = HashMap::new();
    let mut connectivity: HashMap<String, SegmentConnectivity> = HashMap::new();

    // --- PASS 1: Parse and record connectivity + sequences ---
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

            connectivity
                .entry(from.clone())
                .or_default()
                .outgoing
                .insert(to.clone());
            connectivity
                .entry(to.clone())
                .or_default()
                .incoming
                .insert(from.clone());
        }
    }

    // --- PASS 2: Traverse using precomputed in/out counts ---
    let mut current = "s1".to_string();
    let mut visited = HashSet::new();
    let mut unique_path_segments = Vec::new();
    let mut path_count = 1;

    loop {
        if visited.contains(&current) {
            break;
        }
        visited.insert(current.clone());

        let conn = match connectivity.get(&current) {
            Some(c) => c,
            None => break,
        };

        let in_deg = conn.incoming.len();
        let out_deg = conn.outgoing.len();

        if current == "s1" || path_count == 1 {
            unique_path_segments.push(current.clone());
        }

        if out_deg > 1 {
            path_count += out_deg - 1;
        }
        if in_deg > 1 {
            path_count = path_count.saturating_sub(in_deg - 1);
        }

        // Choose lowest-numbered unvisited next segment
        let next = conn
            .outgoing
            .iter()
            .filter(|s| !visited.contains(*s))
            .min_by_key(|s| s[1..].parse::<u64>().unwrap_or(u64::MAX));

        match next {
            Some(n) => current = n.clone(),
            None => break,
        }
    }

    info!("Identified {} unique path segments", unique_path_segments.len());

    // --- Generate kmers ---
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