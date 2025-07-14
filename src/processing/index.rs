use bincode::{DefaultOptions, Options};
use cbl::CBL;
use indicatif::{ProgressBar, ProgressStyle};
use log::*;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use serde::{de::DeserializeOwned, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;
use std::sync::{Arc, Mutex};

use crate::io::utils;
use crate::GeneKmers;

const K: usize = 49;
const PREFIX_BITS: usize = 24;
type T = u128;

fn write_index<S: Serialize, P: AsRef<Path> + Copy>(index: &S, path: P) {
    let output = File::create(path)
        .unwrap_or_else(|_| panic!("Failed to open {}", path.as_ref().to_str().unwrap()));
    let mut writer = BufWriter::new(output);
    info!("Writing the index to {}", path.as_ref().to_str().unwrap());
    DefaultOptions::new()
        .with_varint_encoding()
        .reject_trailing_bytes()
        .serialize_into(&mut writer, &index)
        .unwrap();
}

fn read_index<D: DeserializeOwned, P: AsRef<Path> + Copy>(path: P) -> D {
    let index = File::open(path)
        .unwrap_or_else(|_| panic!("Failed to open {}", path.as_ref().to_str().unwrap()));
    let reader = BufReader::new(index);
    info!(
        "Reading the index stored in {}",
        path.as_ref().to_str().unwrap()
    );
    DefaultOptions::new()
        .with_varint_encoding()
        .reject_trailing_bytes()
        .deserialize_from(reader)
        .unwrap()
}

pub fn build_indexes_for_all_fastas(
    fasta_directory: &Path,
    threads: usize,
    canonical: bool,
    recursive: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    utils::configure_thread_pool(threads);

    let fasta_files =
        utils::find_files_with_extensions(fasta_directory, &["fasta", "fa"], recursive)?;
    let total_files = fasta_files.len();
    if total_files == 0 {
        warn!("No FASTA files found in {:?}", fasta_directory);
        return Ok(());
    }

    info!("Found {} FASTA files to index", total_files);

    // progress bar
    let progress = ProgressBar::new(total_files as u64);
    progress.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%)")
        .unwrap()
        .progress_chars("##-"));

    fasta_files.par_iter().for_each(|fasta_path| {
        info!("Indexing {:?}", fasta_path);

        let result = (|| {
            let mut cbl = if canonical {
                CBL::<K, T>::new_canonical()
            } else {
                CBL::<K, T>::new()
            };

            let mut reader = parse_fastx_file(fasta_path)?;
            while let Some(record) = reader.next() {
                let seqrec = record?;
                cbl.insert_seq(&seqrec.seq());
            }

            let kmers = cbl.count();
            info!(
                "File {:?} contains {} {}{K}-mers",
                fasta_path,
                kmers,
                if canonical { "canonical " } else { "" }
            );

            // Write index next to original file
            let mut index_path = fasta_path.clone();
            index_path.set_extension("cbl");
            write_index(&cbl, &index_path);

            Ok::<_, Box<dyn std::error::Error>>(())
        })();

        if let Err(e) = result {
            warn!("Error indexing {:?}: {}", fasta_path, e);
        }

        progress.inc(1);
    });

    progress.finish_with_message(format!("Indexing complete for all {} files", total_files));
    Ok(())
}

pub fn query_kmers_across_indexes(
    index_directory: &Path,
    mut filtered_kmers: Vec<GeneKmers>,
    threads: usize,
    max_hits: usize,
    recursive: bool,
) -> Result<Vec<GeneKmers>, Box<dyn std::error::Error>> {
    utils::configure_thread_pool(threads);

    let index_files = utils::find_files_with_extensions(index_directory, &["cbl"], recursive)?;
    let total_indexes = index_files.len();
    if total_indexes == 0 {
        warn!("No CBL index files found in {:?}", index_directory);
        return Ok(filtered_kmers);
    }

    info!("Found {} index files to search", total_indexes);

    let mut kmer_to_fk_index: HashMap<String, usize> = HashMap::new();
    let mut kmers: Vec<String> = Vec::new();

    for (i, gene_kmers) in filtered_kmers.iter().enumerate() {
        for probe in &gene_kmers.kmers {
            let kmer = probe.kmer.clone();
            kmer_to_fk_index.insert(kmer.clone(), i);
            kmers.push(kmer);
        }
    }

    info!("Loaded {} kmers from filtered_kmers", kmers.len());

    let results: Arc<Mutex<HashMap<String, Vec<String>>>> = Arc::new(Mutex::new(HashMap::new()));

    let progress = ProgressBar::new(total_indexes as u64);
    progress.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.magenta/blue}] {pos}/{len} ({percent}%)")
        .unwrap()
        .progress_chars("#>-"));

    index_files.par_iter().for_each(|index_path| {
        let result = (|| {
            let mut cbl: CBL<K, T, PREFIX_BITS> = read_index(index_path);
            for kmer in &kmers {
                if cbl.contains_seq(kmer.as_bytes()).iter().any(|&x| x) {
                    let mut res = results.lock().unwrap();
                    res.entry(kmer.clone())
                        .or_default()
                        .push(index_path.to_string_lossy().into_owned());
                }
            }
            Ok::<_, Box<dyn std::error::Error>>(())
        })();

        if let Err(e) = result {
            warn!("Error querying {:?}: {}", index_path, e);
        }

        progress.inc(1);
    });

    progress.finish_with_message("Kmer query complete.");

    let results = results.lock().unwrap();

    for (kmer, files) in results.iter() {
        if let Some(&fk_index) = kmer_to_fk_index.get(kmer) {
            let fk = &mut filtered_kmers[fk_index];
            fk.kmer_hits
                .entry(kmer.clone())
                .or_default()
                .extend(files.clone());
        }
    }

    let filtered = filtered_kmers
        .into_iter()
        .filter(|fk| {
            fk.kmers
                .iter()
                .map(|p| &p.kmer)
                .any(|k| match fk.kmer_hits.get(k) {
                    Some(files) => files.len() <= max_hits,
                    None => true,
                })
        })
        .collect::<Vec<_>>();

    for fk in &filtered {
        for (kmer, files) in fk.kmer_hits.iter().filter(|(_, v)| v.len() <= max_hits) {
            info!(
                "Kmer {} (gene: {}) found in {} index(es):",
                kmer,
                fk.gene,
                files.len()
            );
            for f in files {
                info!("  - {}", f);
            }
        }
    }

    let num_unmatched = kmers.iter().filter(|k| !results.contains_key(*k)).count();
    info!(
        "{} of {} kmers had no hits in any index.",
        num_unmatched,
        kmers.len()
    );

    if filtered.is_empty() {
        warn!("All kmers were filtered out - no kmers matched the criteria");
    }

    Ok(filtered)
}
