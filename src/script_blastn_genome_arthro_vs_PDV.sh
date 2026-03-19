#!/usr/bin/env bash

##############################################################################
# CONFIGURATION — edit these variables before running
##############################################################################

GENOME_IDS="id_genome_protostomia.txt"
PDV_FASTA="pdv.fasta"
N_SUBSETS=10
DIR_RAW="raw/"
DIR_PROCESSED="processed/"
DIR_LOGS="logs/"
THREADS=8

##############################################################################
# SETUP
##############################################################################

mkdir -p "$DIR_RAW" "$DIR_PROCESSED" "$DIR_LOGS"

# Generate zero-padded subset IDs: 00, 01, ..., N_SUBSETS-1
SUBSET_IDS=()
for i in $(seq 0 $(( N_SUBSETS - 1 ))); do
    SUBSET_IDS+=( $(printf '%02d' $i) )
done

##############################################################################
# STAGE 1 — Split id_genome_protostomia.txt into N_SUBSETS subsets
##############################################################################

echo "Splitting genome ID list into $N_SUBSETS subsets..."

total=$(wc -l < "$GENOME_IDS")
lines_per_chunk=$(( (total + N_SUBSETS - 1) / N_SUBSETS ))

split -l "$lines_per_chunk" -d -a 2 --additional-suffix='.txt' \
    "$GENOME_IDS" "${DIR_RAW}subset-" 2> "${DIR_LOGS}split_genome_ids.log"

# Ensure all expected subset files exist (create empty ones if fewer chunks)
for num in "${SUBSET_IDS[@]}"; do
    src="${DIR_RAW}subset-${num}.txt"
    if [ ! -f "$src" ]; then
        touch "$src"
    fi
done

echo "Done."

##############################################################################
# STAGE 2 — BLAST each subset genome against PDV
##############################################################################

echo "Running BLAST for each subset..."

for num in "${SUBSET_IDS[@]}"; do
    subset="${DIR_RAW}subset-${num}.txt"
    blast_out="${DIR_PROCESSED}ALL_outputblastn_${num}.txt"
    url_log="${DIR_LOGS}url_genome_${num}.txt"
    log="${DIR_LOGS}blast_${num}.log"

    echo "  Processing subset $num..."

    tmpdir=$(mktemp -d)

    (
        cd "$tmpdir"

        while read -r acc; do
            echo "NEW GENOME: $acc"

            if [[ $acc = GCF* ]]; then
                datatype=FtpPath_RefSeq
            else
                datatype=FtpPath_GenBank
            fi
            echo "$datatype"

            base_acc=$(echo "$acc" | awk -F "." '{print $1}')
            echo "$base_acc"

            url=$(esearch -db assembly -query "$acc" </dev/null \
                  | esummary \
                  | xtract -pattern DocumentSummary -element "$datatype" \
                  | grep "$base_acc" \
                  | awk '{sub("ftp://", ""); print}')

            dir=$(echo "$url" | grep -o 'GC._.*')
            fname="${url}/${dir}_genomic.fna.gz"
            echo "$fname" >> "$url_log"
            echo "Downloading $acc from $fname"

            wget "$fname" --quiet
            gunzip *.gz* 2>/dev/null || true

            makeblastdb -in *.fna -dbtype nucl
            blastn -task blastn \
                   -query "$PDV_FASTA" \
                   -db *.fna \
                   -outfmt '6 qaccver qlen saccver stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
                   -num_threads "$THREADS" \
                   -evalue 0.0001 \
                   -out "outputblastn_${base_acc}"

            rm *.fna*

            find . -iname "outputblastn_${base_acc}" -print0 | xargs -0 awk '
                FNR==1{ sub(/\.txt$/, "", FILENAME) }
                FNR==1{ sub("./",    "", FILENAME) }
                FNR==1{ sub("outputblastn_", "", FILENAME) }
                { print $0 "\t" FILENAME }' \
                >> "$blast_out"

            rm "outputblastn_${base_acc}"

        done < "$subset"

    ) 2>> "$log"

    rm -rf "$tmpdir"
    touch "$blast_out"   # avoid missing-file errors if no hits

    echo "  Subset $num done."
done

echo "Done."

##############################################################################
# STAGE 3 — Concatenate all per-subset BLAST results
##############################################################################

echo "Concatenating all BLAST results..."

final_output="${DIR_RAW}outputblastn_arthro_vs_PDV_ALL_11.txt"
log="${DIR_LOGS}concatenate_blast.log"

cat "${DIR_PROCESSED}"ALL_outputblastn_*.txt > "$final_output" 2> "$log"

echo "Done."
echo "Final output: $final_output"
