use log::*;
use visiogen::FilteredKmers;

use crate::cli::GraphArgs;
use crate::kmer;

use std::collections::HashMap;
use std::io::BufRead;

pub struct Gfa {
    pub segments: Vec<Segment>,
    pub links: Vec<Link>,
    pub paths: Vec<GfaPath>,
}

impl Gfa {
    /// Return segment names that appear exactly once in all paths (core)
    pub fn core_segments(&self) -> Vec<String> {
        let path_count = self.paths.len();
        let mut segment_in_path_counts: HashMap<String, usize> = HashMap::new();

        for path in &self.paths {
            let mut segment_seen: HashMap<String, usize> = HashMap::new();

            for (segment_name, _) in &path.segments {
                *segment_seen.entry(segment_name.clone()).or_insert(0) += 1;
            }

            for (seg, count) in segment_seen {
                if count == 1 {
                    *segment_in_path_counts.entry(seg).or_insert(0) += 1;
                }
            }
        }

        segment_in_path_counts
            .into_iter()
            .filter_map(
                |(seg, count)| {
                    if count == path_count {
                        Some(seg)
                    } else {
                        None
                    }
                },
            )
            .collect()
    }

    /// Return full Segment structs instead of just names
    pub fn core_segment_structs(&self) -> Vec<&Segment> {
        let core_names = self.core_segments();
        let name_set: std::collections::HashSet<_> = core_names.iter().collect();

        self.segments
            .iter()
            .filter(|seg| name_set.contains(&seg.name))
            .collect()
    }
}

enum GfaLine {
    Segment(Segment),
    Link(Link),
    Path(GfaPath),
    Other(String),
}

#[derive(Debug)]
pub struct Segment {
    pub name: String,
    pub sequence: String,
}

#[derive(Debug)]
pub struct Link {
    from: String,
    from_orient: char,
    to: String,
    to_orient: char,
    overlap: String,
}

#[derive(Debug)]
pub struct GfaPath {
    name: String,
    segments: Vec<(String, char)>,
    overlaps: Vec<String>,
}

fn parse_line(line: &str) -> Option<GfaLine> {
    let fields: Vec<&str> = line.split('\t').collect();

    match fields.get(0)? {
        &"S" => Some(GfaLine::Segment(Segment {
            name: fields.get(1)?.to_string(),
            sequence: fields.get(2)?.to_string(),
        })),
        &"L" => Some(GfaLine::Link(Link {
            from: fields.get(1)?.to_string(),
            from_orient: fields.get(2)?.chars().next()?,
            to: fields.get(3)?.to_string(),
            to_orient: fields.get(4)?.chars().next()?,
            overlap: fields.get(5)?.to_string(),
        })),
        &"P" => {
            let name = fields.get(1)?.to_string();
            let segments: Vec<(String, char)> = fields
                .get(2)?
                .split(',')
                .map(|s| {
                    let (seg, orient) = s.split_at(s.len() - 1);
                    (seg.to_string(), orient.chars().next().unwrap())
                })
                .collect();

            let overlaps = fields.get(3)?.split(',').map(|o| o.to_string()).collect();

            Some(GfaLine::Path(GfaPath {
                name,
                segments,
                overlaps,
            }))
        }
        _ => Some(GfaLine::Other(line.to_string())),
    }
}

pub fn parse_gfa_file(path: &str) -> std::io::Result<Gfa> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);

    let mut segments = Vec::new();
    let mut links = Vec::new();
    let mut paths = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        if let Some(parsed) = parse_line(&line) {
            match parsed {
                GfaLine::Segment(s) => segments.push(s),
                GfaLine::Link(l) => links.push(l),
                GfaLine::Path(p) => paths.push(p),
                GfaLine::Other(_) => (),
            }
        }
    }

    Ok(Gfa {
        segments,
        links,
        paths,
    })
}

pub fn run_graph_mode(graph_args: &GraphArgs, kmer_size: usize) -> Vec<FilteredKmers> {
    let graph = parse_gfa_file(&graph_args.gfa_path).expect("Failed to read GFA file");

    let filtered_kmers: Vec<FilteredKmers> = graph
        .core_segment_structs()
        .iter()
        .map(|segment| FilteredKmers {
            gene: segment.name.clone(),
            start: 1,
            end: 1 + segment.sequence.len() as u64,
            kmers: kmer::tile_segment(&segment.sequence, 1 as usize, kmer_size),
            strand: "+".to_string(),
            kmer_hits: HashMap::new(),
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