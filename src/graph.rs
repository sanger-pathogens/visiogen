use log::*;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use visiogen::FilteredKmers;
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::Direction;

use crate::cli::GraphArgs;
use crate::kmer;

#[derive(Debug)]
struct Segment {
    id: String,
    sequence: String,
    start: u64,
}

#[derive(Debug)]
struct EdgeAttributes {
    sr: u32,
}

fn parse_segments(path: &str) -> (HashMap<String, Segment>, Vec<String>) {
    let file = File::open(path).expect("Failed to open GFA file");
    let reader = BufReader::new(file);

    let mut segments: HashMap<String, Segment> = HashMap::new();
    let mut segment_order: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if !line.starts_with('S') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        let seg_id = parts[1].to_string();
        let seq = parts[2].to_string();
        let mut so = 0u64;

        for field in &parts[3..] {
            if field.starts_with("SO:i:") {
                so = field[5..].parse().unwrap_or(0);
            }
        }

        segments.insert(
            seg_id.clone(),
            Segment {
                id: seg_id.clone(),
                sequence: seq,
                start: so,
            },
        );
        segment_order.push(seg_id);
    }

    (segments, segment_order)
}

fn parse_links(path: &str) -> Vec<(String, String, u32)> {
    let file = File::open(path).expect("Failed to open GFA file");
    let reader = BufReader::new(file);
    let mut links = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if !line.starts_with('L') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 6 {
            continue;
        }

        let from = parts[1].to_string();
        let to = parts[3].to_string();

        let mut sr_val = 0u32;
        for field in &parts[5..] {
            if field.starts_with("SR:i:") {
                sr_val = field[5..].parse().unwrap_or(0);
                break;
            }
        }

        links.push((from, to, sr_val));
    }

    links
}

fn build_graph(
    segments: &HashMap<String, Segment>,
    links: &[(String, String, u32)],
) -> (DiGraph<String, EdgeAttributes>, HashMap<String, NodeIndex>) {
    let mut graph = DiGraph::<String, EdgeAttributes>::new();
    let mut node_map = HashMap::new();

    // Add nodes
    for seg_id in segments.keys() {
        let node = graph.add_node(seg_id.clone());
        node_map.insert(seg_id.clone(), node);
    }

    // Add edges
    for (from, to, sr) in links {
        if let (Some(&from_node), Some(&to_node)) = (node_map.get(from), node_map.get(to)) {
            graph.add_edge(from_node, to_node, EdgeAttributes { sr: *sr });
        }
    }

    (graph, node_map)
}

fn traverse_bubble_depth(
    graph: &DiGraph<String, EdgeAttributes>,
    node_map: &HashMap<String, NodeIndex>,
    start_id: &str,
) -> Vec<String> {
    let mut segments = Vec::new();
    let mut visited = HashSet::new();
    let mut current_node = match node_map.get(start_id) {
        Some(&n) => n,
        None => return segments,
    };

    let mut bubble_depth = 0;
    segments.push(graph[current_node].clone());
    visited.insert(current_node);

    // Track which edges we've actually traversed
    let mut traversed_edges = HashSet::new();

    loop {
        // Get all outgoing edges
        let outgoing_edges: Vec<_> = graph.edges_directed(current_node, Direction::Outgoing).collect();
        
        // Partition into SR0 and non-SR0 edges
        let (sr0_edges, non_sr0_edges): (Vec<_>, Vec<_>) = 
            outgoing_edges.into_iter().partition(|e| e.weight().sr == 0);

        // Find the next SR0 edge to follow (should be only one unvisited)
        let next_sr0_edge = sr0_edges.iter()
            .find(|e| !visited.contains(&e.target()));

        // Count potential bubble branches (non-SR0 edges)
        let new_bubbles = non_sr0_edges.iter()
            .filter(|e| !visited.contains(&e.target()))
            .count();

        // Update bubble depth when we have branching paths
        if new_bubbles > 0 {
            bubble_depth += new_bubbles;
            debug!("Entering bubble at {}: depth increased to {}", 
                  graph[current_node], bubble_depth);
        }

        match next_sr0_edge {
            Some(edge) => {
                let next_node = edge.target();
                traversed_edges.insert((current_node, next_node));
                visited.insert(next_node);
                
                // Check for any incoming edges that represent bubble merges
                let incoming_non_sr0: Vec<_> = graph.edges_directed(next_node, Direction::Incoming)
                    .filter(|e| e.weight().sr != 0)
                    .collect();

                if !incoming_non_sr0.is_empty() {
                    // When we're in a complex bubble (depth >= 2), investigate properly
                    let merges_to_subtract = if bubble_depth >= 2 {
                        investigate_bubble(graph, next_node, &traversed_edges)
                    } else {
                        incoming_non_sr0.len()
                    };

                    bubble_depth = bubble_depth.saturating_sub(merges_to_subtract);
                    debug!("Merging bubble at {}: subtracted {}, depth now {}", 
                          graph[next_node], merges_to_subtract, bubble_depth);
                }

                // Only add to segments if we're not in a bubble
                if bubble_depth == 0 {
                    segments.push(graph[next_node].clone());
                    debug!("Adding segment: {}", graph[next_node]);
                } else {
                    debug!("Skipping segment {} (depth {})", 
                          graph[next_node], bubble_depth);
                }
                
                current_node = next_node;
            },
            None => {
                // No SR0 edges available - try non-SR0 edges if we're in a bubble
                if bubble_depth > 0 {
                    if let Some(edge) = non_sr0_edges.iter()
                        .find(|e| !visited.contains(&e.target())) 
                    {
                        let next_node = edge.target();
                        info!("Following bubble path to {}", graph[next_node]);
                        traversed_edges.insert((current_node, next_node));
                        visited.insert(next_node);
                        current_node = next_node;
                    } else {
                        break;
                    }
                } else {
                    break; // No SR0 edges and not in a bubble
                }
            }
        }
    }

    segments
}

fn investigate_bubble(
    graph: &DiGraph<String, EdgeAttributes>,
    merge_node: NodeIndex,
    traversed_edges: &HashSet<(NodeIndex, NodeIndex)>,
) -> usize {
    let mut visited_nodes = HashSet::new();
    let mut visited_edges = HashSet::new();
    let mut queue = VecDeque::new();
    let mut merge_count = 0;

    // Start with all incoming non-SR0 edges to the merge node
    for edge in graph.edges_directed(merge_node, Direction::Incoming) {
        if edge.weight().sr != 0 {
            queue.push_back((edge.source(), edge.id()));
        }
    }

    while let Some((node, edge_id)) = queue.pop_front() {
        if visited_nodes.contains(&node) {
            continue;
        }
        visited_nodes.insert(node);
        visited_edges.insert(edge_id);

        // Check ALL incoming edges to this node
        for edge in graph.edges_directed(node, Direction::Incoming) {
            if edge.weight().sr == 0 {
                // Found connection to main SR0 path - count this edge
                if traversed_edges.contains(&(edge.source(), node)) {
                    merge_count += 1;
                }
            } else if !visited_edges.contains(&edge.id()) {
                // Continue tracing back through non-SR0 edges we haven't processed
                queue.push_back((edge.source(), edge.id()));
            }
        }
    }

    debug!("Bubble investigation at {}: found {} merging edges", 
          graph[merge_node], merge_count);
    merge_count
}


pub fn run_graph_mode(graph_args: &GraphArgs, kmer_size: usize) -> Vec<FilteredKmers> {
    let (segments, _segment_order) = parse_segments(&graph_args.gfa_path);
    let links = parse_links(&graph_args.gfa_path);

    let (graph, node_map) = build_graph(&segments, &links);

    let chosen_segments = traverse_bubble_depth(&graph, &node_map, "s1");

    let filtered_kmers: Vec<FilteredKmers> = chosen_segments
        .iter()
        .filter_map(|seg_id| {
            let segment = segments.get(seg_id)?;
            Some(FilteredKmers {
                gene: seg_id.clone(),
                start: segment.start,
                end: segment.start + segment.sequence.len() as u64,
                kmers: kmer::tile_segment(&segment.sequence, segment.start as usize, kmer_size),
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
