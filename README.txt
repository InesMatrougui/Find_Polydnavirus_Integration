# PDV Integration in Host Genome Analysis
# PDV

## Overview

This R script analyzes polydnavirus (PDV) integration patterns in protostome genomes using BLAST search results. Protostome genomes are available on NCBI: https://www.ncbi.nlm.nih.gov/datasets/genome/. The analysis focuses on identifying genomic integrations mediated by Host Integration Motifs (HIMs) while filtering out false positives.

This project includes an R script that processes data using a combination of:
- R code (organized into executable chunks)
- External command-line tools (Linux-based)

**System Requirements:** Linux or Unix-based OS

**Important:** The workflow relies on running each R chunk individually, with manual execution of external commands in between chunks as specified in the script.

## Author
Ines Matrougui  
Date: 2025-09-04

## Purpose

To identify and characterize polydnavirus circle integrations in host genomes by:
- Filtering BLAST hits to retain those falling within Host Integration Motifs (HIMs)
- Removing false positives (low complexity regions, contamination)
- Analyzing taxonomic distribution of integrations
- Identifying orthologous integrations across species
- Performing phylogenetic analysis of integration patterns

## Input Data

### Required Files
Place all input files in the `data/raw/` directory:

- `outputblastn_arthro_vs_PDV_ALL_11.txt` - Tabular BLASTN output from genome searches
- `HIMcoord.txt` - Coordinates of HIM boundaries in PDV segments
- `seg_new_HIM.txt` - Coordinates of newly identified HIM boundaries
- `Segments.txt` - Names of PDV segments (no header)
- `PDV_wasp_sp.txt` - PDV segments with corresponding wasp species
- `id_genome_protostomia.txt` - Protostome genome identifiers from NCBI
- `taxonomy_genome.txt` - Taxonomic information (generated via efetch)
- `coor_LCR_in_hit_seg_06_12.txt` - Low complexity region coordinates
- `toMerge_merged_expand_all.txt` - Merged coordinates from bedtools
- `position_hit_merged.txt` - Position data after overlap removal
- `taxResult05_flank_tax_per_contig_read.tsv` - Metaeuk taxonomic assignment results
- `outputmegablast_ortho_flank.txt` - Megablast results for orthology analysis

## Dependencies

### R Packages
```r
library(data.table)
library(dplyr) 
library(tidyverse)
library(ggplot2)
library(igraph)
library(callr)
```

### External Tools Required
- **BLASTN** - For sequence similarity searches
- **bedtools** - For coordinate merging and intersection operations
- **RepeatMasker** - For low complexity region identification  
- **seqtk** - For sequence extraction
- **efetch/esearch** - For sequence and taxonomy retrieval
- **metaeuk** - For taxonomic assignment of contigs
- **mmseqs2** - For clustering and database creation
- **iqtree** - For phylogenetic analysis

## Directory Structure

```
project/
├── data/
│   ├── raw/           # Input files
│   └── processed/     # Intermediate processed files
├── output/            # Final results and tables
├── figures/           # Generated plots and phylogenies
└── src/               # Additional R scripts for phylogeny
```

## Workflow Overview

### 1. Data Import and Preprocessing
- Load BLASTN results and filter by e-value (< 0.0001) and alignment length (> 164 bp)
- Import HIM coordinates and PDV segment information
- Clean column names and add indexing

### 2. HIM-Mediated Integration Analysis
- Filter hits involving segments with known HIMs
- Identify hits spanning HIM boundaries (unexpected if integration is HIM-mediated)
- Categorize hits by junction type (J1 or J2)
- Process new HIM data and merge with existing results

### 3. Relaxed Stringency Analysis
- Accept hits within 100bp of HIM boundaries
- Classify hits as inside/outside HIM regions
- Combine strict and relaxed criteria results

### 4. Taxonomic Classification
- Extract species names from sequence descriptions
- Link genome accessions to taxonomic information
- Filter out parasitoid wasps (Braconidae, Ichneumonidae)

### 5. False Positive Removal
- Remove hits overlapping with low complexity regions (LCRs)
- Merge overlapping coordinates to remove redundancy
- Filter based on taxonomic assignment of flanking regions (retain only eukaryotic)

### 6. Integration Validation
- Extract flanking sequences for taxonomic validation
- Use metaeuk for taxonomic assignment of contigs
- Retain only integrations in eukaryotic contigs

### 7. Phylogenetic Analysis
- Select representative species for phylogeny
- Analyze integration patterns across taxa
- Create species phylogenies showing integration distribution
- Perform junction-specific (J1/J2) phylogenetic analysis

### 8. Orthology Analysis
- Identify similar integrations across species
- Use sequence similarity and genomic context
- Build networks of orthologous integrations
