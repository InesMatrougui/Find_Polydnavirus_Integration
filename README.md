# PDV Integration Analysis in Protostome Genomes

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
**Please note that the `script_blastn_genome_arthro_vs_PDV.sh` can generate the `outputblastn_arthro_vs_PDV_ALL_11.txt`**
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

### Steps to take beforehand `[Bash]`
- Split genome IDs 
Splits `id_genome_protostomia.txt` into 8 equal subsets (subset-00.txt … subset-07.txt) using split. This enables parallel BLAST jobs.
- BLAST subsets against PDV
For each genome accession in a subset:
    Determines whether to use the RefSeq or GenBank FTP path based on the accession prefix (GCF_ → RefSeq, GCA_ → GenBank).
    Downloads the genome FASTA via esearch/esummary/xtract + wget.
    Builds a BLAST nucleotide database with makeblastdb.
    Runs blastn with the PDV HIM sequences as query (-evalue 0.0001, 8 threads).
    Appends the genome accession ID as a final column to each hit line.
    Deletes the genome FASTA to save disk space.
All 8 subset jobs run in parallel (one Snakemake job per subset).
- Concatenate BLAST results
Concatenates the 8 per-subset BLAST output files into `outputblastn_arthro_vs_PDV_ALL_11.txt`, which is the input expected by the R analysis.

### 1. Data Import and Preprocessing `[R]`
- Load BLASTN results and filter by e-value (< 0.0001) and alignment length (> 164 bp)
- Import HIM coordinates and PDV segment information
- Clean column names and add indexing

### 2. HIM-Mediated Integration Analysis `[R]`
- Filter hits involving segments with known HIMs
- Identify hits spanning HIM boundaries (unexpected if integration is HIM-mediated)
- Categorize hits by junction type (J1 or J2)
- Process new HIM data and merge with existing results

### 3. Fetch contig sequences `[Bash]`
Downloads all contig FASTA sequences from NCBI nucleotide using esearch + efetch, producing contig.fasta.

### 4. Relaxed Stringency Analysis `[R]`
- Accept hits within 100bp of HIM boundaries
- Classify hits as inside/outside HIM regions
- Combine strict and relaxed criteria results

### 5. Taxonomic Classification `[R]`
- Extract species names from sequence descriptions
- Link genome accessions to taxonomic information
  #### 5.a Fetch taxonomy `[Bash]`
Retrieves full taxonomic lineages for all TaxIDs using `efetch -format xml` and `xtract`, producing `taxonomy_genome.txt` used to annotate hits with order, family, and genus.
- Filter out parasitoid wasps (Braconidae, Ichneumonidae)

### 6 False Positive Removal `[R]`
- Remove hits overlapping with low complexity regions (LCRs)
  #### 6.a False Positive RemovalStage with RepeatMasker `[Bash]`
Annotates low complexity regions (LCRs) in the PDV query sequences using RepeatMasker (-low flag restricts to LCR and simple repeats only). These represent false-positive hits arising from sequence composition rather than true homology.
Identifies the coordinates of LCRs that fall within BLAST hits with `bedtools intersect`, producing coor_LCR_in_hit_seg_06_12.txt.

- Merge overlapping coordinates to remove redundancy
  #### 6.b bedtools merge `[Bash]`
Merges overlapping hit coordinates strand-specifically to remove redundancy caused by multiple similar PDV segments (e.g. Hd12 and Hd16, or BV10 and BV17) aligning to the same host locus. dos2unix is applied first to prevent carriage-return issues.

- Filter based on taxonomic assignment of flanking regions (retain only eukaryotic)
Extracts 200,000 bp flanking windows around each integration for metaeuk taxonomic assignment.
  #### 6.c — Flank extraction and metaeuk `[Bash]`
1. `bedtools merge` collapses overlapping flank windows.
2. `seqtk subseq` extracts flank sequences from `contig.fasta`.
3. `metaeuk easy-predict` predicts proteins in flank sequences against the NR database.
4. `metaeuk taxtocontig` assigns taxonomy to each contig based on the taxonomy of its predicted proteins (majority vote, LCA mode 4).

### 7. Integration Validation `[R]`
- Retain only integrations in eukaryotic contigs
- Removes alternate haplotype assemblies (keeps one assembly per species).
- Computes summary statistics (number of integrations, species, segments, identity, length).
- Standardises PDV circle names across homologous segments (e.g. Hd27 ↔ KF156228.1_DsIV_15 → IV27).
- Writes the final integration table, species phylogeny list, integration counts, and position files for phylogeny and orthology.
  
### 8. Phylogenetic Analysis `[R]`
- Select representative species for phylogeny
- Analyze integration patterns across taxa
- Create species phylogenies showing integration distribution
- Perform junction-specific (J1/J2) phylogenetic analysis
1. `seqtk subseq` extracts Hd27/Ds15 J1 and J2 sequences for phylogenetic analysis.
2. **Alignment (muscle/Geneious) is performed externally** 
3. `trimal -automated1` trims the alignment.
4. `iqtree -m MFP -bb 1000` infers the maximum-likelihood phylogeny.

### 9. Orthology Analysis `[R]`
- Identify similar integrations across species
1. `bedtools merge` collapses ortho flank windows.
2. `seqtk subseq` extracts ortho flank sequences.
3. `blastn -task megablast` searches flanks against themselves to find orthologous integrations.
- Imports megablast results and annotates both the query and subject sides with integration coordinates and species metadata.
- Filters for inter-species hits where the alignment covers the integration region and extends meaningfully into the flanking contig (>50% flank, >90% hit coverage).
- Exports orthologous integration coordinates and writes the network graph via `igraph`.
- Calls `phylogeny_species.R` via `callr` to generate SVG figures.

---

## Notes
- the script 
