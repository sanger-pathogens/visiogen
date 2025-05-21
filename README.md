# Visiogen

![Rust Version](https://img.shields.io/badge/Rust-nightly--2024--05--31-blue?style=flat-square)
![License](https://img.shields.io/badge/License-AGPL--3.0-blue?style=flat-square)
[![Issues](https://img.shields.io/github/issues/sanger-pathogens/visiogen)](https://github.com/sanger-pathogens/visiogen/issues)

**Visiogen** is a high-performance tool for producing probes using **kmers** from genomic annotations/graphs. It supports GFF-based gene input, graph-based input, and kmer indexing with customizable filtering using the CBL libary.

Check out [CBL](https://github.com/imartayan/CBL/tree/main) for more details.

---

## 🔧 Features

* Generate probes from GFA graphs or GFF files
* Filter kmers based on GC content, center base, or other sequence characteristics
* Build off target kmer indexes from FASTA inputs
* Query kmers against index files to assess off-target risk

---

## 🚀 Installation

You’ll need:

* Rust toolchain (`nightly-2024-05-31`)
* System dependencies:

```bash
sudo apt install -y libstdc++-12-dev libclang-dev
```

then after cloning the repo you can build with

```bash
RUN cargo +nightly-2024-05-31 build --release
```

### 🐳 Docker

Two Dockerfiles are provided:

* **Ubuntu 22.04**: `Dockerfile`
* **Ubuntu 20.04**: `Dockerfile-ubuntu-20`

Build example:

```bash
docker build -f Dockerfile -t visiogen .
```
---

## 🧬 Example Usage

### 🔹 GFF Mode: Extract kmers from specific genes

```bash
visiogen gff -f input.fa -a annotation.gff -g gene1,gene2,gene3
```

### 🔹 Graph Mode: Generate probes from a GFA assembly graph

```bash
visiogen graph -g graph.gfa -t 0.95
```

### 🔹 Build Mode: Index FASTA files for off-target querying

```bash
visiogen build -i fasta_dir
```

### 🔹 After indexing the fasta_dir you can use it as an off target database

```bash
visiogen graph -g graph.gfa -t 0.95 -i fasta_dir
```

---

## ⚙️ Global Arguments

| Flag                         | Description                                            |
| ---------------------------- | ------------------------------------------------------ |
| `-t, --threads`              | Number of threads to use (default: all cores)          |
| `-i, --off_target_directory` | Directory of FASTA/index files to scan for off-targets |
| `--max_hits`                 | Max index hits per kmer to retain (default: 5)         |
| `-r, --recursive`            | Recursively scan directories for index files           |

---

### Kmer Options

| Flag                | Description                                       |
| ------------------- | ------------------------------------------------- |
| `-k, --kmer_size`   | Length of kmers (default: 50)                     |
| `-b, --center_base` | Center base (e.g., G) to constrain selection      |
| `-l, --min_gc`      | Minimum GC content (default: 44)                  |
| `-m, --max_gc`      | Maximum GC content (default: 72)                  |
| `--allow_outside`   | Allow kmers outside target genes (default: false) |
| `--skip_gc`         | Disable GC filtering                              |

---

## 📦 Subcommands

### `gff`

Extract kmers from a genome FASTA using GFF annotation.

Required:

* `-f <FASTA>`: Genome sequence
* `-a <GFF>`: Gene annotation
* `-g <genes>`: Comma-separated list of gene IDs

### `build`

Create `.cbl` kmer index files from a directory of FASTA files.

Optional:

* `-c, --canonical`: Use canonical kmers (default: true)

### `graph`

Generate kmers from a GFA-format assembly graph.

* `-g <GFA>`: Path to `.gfa` graph
* `-t <threshold>`: Core segment threshold (default: 0.95)

---

## 📤 Output

### Each kmer is tracked with:

* Associated gene, strand, and coordinates
* Presence/absence in index files (off-targets)
* Filtering status based on GC content, hits, etc.

### Logging reports:

* Which kmers were found in which indexes
* Which were discarded due to too many hits
* Which were missing from all indexes

---

## 🧼 Logging Examples

```
13:31:06 [INFO] Found 4 strains; retaining segments with SR:i ≥ 4 (threshold = 0.95)
13:31:06 [INFO] Total segments available: 1008
13:31:06 [INFO] Segments passing strain threshold: 93
13:31:06 [INFO] Generated kmers for 93 segments (total raw kmers: 420928, avg per segment: 4526.11)
13:31:06 [INFO] Using default number of threads (all logical CPUs)
13:31:06 [INFO] Found 12 index files to search
13:31:06 [INFO] Loaded 62264 kmers from filtered_kmers
13:31:06 [INFO] Reading the index stored in fastas/GCA_037443265.cbl
13:31:07 [INFO] Kmer CCGCCGACTGCCCCAAATTCACGCTGTCAGAGACGGTATTGGCGATATCG (gene: s984) found in 1 index(es):
13:31:07 [INFO]   - fastas/GCA_000069185.cbl
```

---

## 🛣️ Roadmap

Planned features and improvements:

### Graph Mode Enhancements
Stricter support thresholds: Link segment inclusion thresholds to actual strain support.
Bubble resolution: Interpret graph bubbles as SNPs or structural variation to increase number of usable segments.
Annotation: Support importing annotations from a reference genome and adding to the graph.

### 🧹 Logging & UX
Improved logging: Reduce verbosity and group messages.

### 🗃️ Indexing Improvements
Single-file index format: Consolidate multiple .cbl files into a unified index.

### 🦀 Rust & Build System
Migrate from nightly toolchain to stable Rust once key dependencies allow/are changed.

Offer precompiled binaries or automated builds for common platforms (e.g., Linux, macOS).

## 📫 Feedback / Contributions

Issues welcome! This was a learning project for me and any feedback or ideas are welcome 😊

> **⚠️ WARNING:** Always manually BLAST the selected probes before ordering them to verify specificity and avoid off-target binding.
