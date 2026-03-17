# ==============================================================
# FUNCTION TO ORDER DNA SEQUENCES BASED ON STRAND ORIENTATION
# ==============================================================

library(data.table)
library(dplyr)
library(Biostrings)
library(tidyr)

# ----------------------
# Define Input/Output Paths
# ----------------------
dir_In_processed <- "data/processed/"
dir_Output <- "output/"

# ----------------------
# Define the function
# ----------------------
process_fasta_by_strand <- function(fasta_path, position_df, output_fasta_path) {
  
  # Step 1: Load DNA sequences from FASTA
  dna <- readDNAStringSet(paste0(dir_In_processed,fasta_path))
  
  # Step 2: Extract and clean sequence metadata
  seq_name <- names(dna)
  ref_seq <- sub(" .*", "", seq_name)  # Remove any extra info after first space
  sequence <- paste(dna)
  df <- data.frame(ref_seq, seq_name, sequence)
  
  # Step 3: Format the position table (standardize column names)
  position_df <- read.table(paste0(dir_Output, position_df), header = FALSE)
  setnames(position_df, c("ref", "start", "end", "stran"))
  
  # Step 4: Format coordinates to match FASTA IDs
  id <- position_df %>%
    mutate(start = start + 1) %>%
    unite("Ref", ref:start, sep = ":") %>%
    unite("ref_seq", Ref:end, sep = "-") %>%
    select(ref_seq, stran)
  
  # Step 5: Split by strand
  id_minus <- filter(id, stran == '-')
  id_plus  <- filter(id, stran == '+')
  
  # Step 6: Process minus strand (reverse complement)
  minus <- merge(id_minus, df, by = "ref_seq")
  filtered_ref <- dna[names(dna) %in% minus$seq_name]
  dna_minus <- reverseComplement(filtered_ref)
  
  # Step 7: Process plus strand (keep as-is)
  plus <- merge(df, id_plus, by = "ref_seq")
  dna_plus <- dna[names(dna) %in% plus$seq_name]
  
  # Step 8: Combine and write sequences
  ordered_dna <- c(dna_plus, dna_minus)
  writeXStringSet(ordered_dna, paste0(dir_Output, output_fasta_path))
  
}

# ==============================================================

# ----------------------
# For J1 (phylo J1/J2)
# ----------------------
process_fasta_by_strand(
  fasta_path = "hit_Ds15_Hd27_no_Lepido_J1.fasta",
  position_df = "position_hit_Ds15_Hd27_no_Lepido_J1.txt",
  output_fasta_path = "hit_Ds15_Hd27_no_Lepido_J1_ordered.fasta"
)

# ----------------------
# For J2 (phylo J1/J2)
# ----------------------
process_fasta_by_strand(
  fasta_path = "hit_Ds15_Hd27_no_Lepido_J2.fasta",
  position_df = "position_hit_Ds15_Hd27_no_Lepido_J2.txt",
  output_fasta_path = "hit_Ds15_Hd27_no_Lepido_J2_ordered.fasta"
)


# ==============================================================

# ----------------------
# For J1 (clustering)
# ----------------------
process_fasta_by_strand(
  fasta_path = "hit_CLUS_J1.fasta",
  position_df = "position_hit_CLUS_J1.txt",
  output_fasta_path = "hit_CLUS_J1_ordered.fasta"
)

# ----------------------
# For J2 (clustering)
# ----------------------
process_fasta_by_strand(
  fasta_path = "hit_CLUS_J2.fasta",
  position_df = "position_hit_CLUS_J2.txt",
  output_fasta_path = "hit_CLUS_J2_ordered.fasta"
)

# ----------------------
# For J1 (clustering) IV27
# ----------------------
process_fasta_by_strand(
  fasta_path = "hit_CLUS_IV27_J1.fasta",
  position_df = "position_hit_CLUS_IV27_J1.txt",
  output_fasta_path = "hit_CLUS_IV27_J1_ordered.fasta"
)

# ----------------------
# For J2 (clustering) IV27
# ----------------------
process_fasta_by_strand(
  fasta_path = "hit_CLUS_IV27_J2.fasta",
  position_df = "position_hit_CLUS_IV27_J2.txt",
  output_fasta_path = "hit_CLUS_IV27_J2_ordered.fasta"
)



# removing duplicate : seqkit rmdup HIM_arthro_ordered.fasta > HIM_arthro_ordered_uniq.fasta
# Clustering sequences : mmseqs easy-cluster /mnt/65To/Ines/PDV_arthro/HIM_arthro_ordered.fasta clusterRes tmp --min-seq-id 0.9 -c 0.2 --cov-mode 1