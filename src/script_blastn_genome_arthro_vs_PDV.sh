num=aa
acc_file=subset-${num}

## For each ID from acc_file dowload the genome and blast against PDV
cat ${acc_file} | while read -r acc ; do
    ## To decide Refseq or Genebank
    echo "NEW GENOME"
    if [[ $acc = GCF* ]];then datatype=FtpPath_RefSeq;else datatype=FtpPath_GenBank;fi
    echo ${datatype}
    ## base_acc = GenBank ID without .end
    base_acc=`echo ${acc}|awk -F "." '{print $1}'`
    echo ${base_acc}
    ## get url Genbank
    url=`esearch -db assembly -query ${acc} </dev/null|esummary|xtract -pattern DocumentSummary -element ${datatype}|grep "$base_acc"|awk '{sub("ftp://", ""); print}'`
    dir=`echo ${url}|grep -o 'GC._.*'`
    fname=${url}/${dir}_genomic.fna.gz
    echo "$fname" >> url_genome_${num}.txt
    echo "Downloading $acc from $fname"
    ## download genome
    wget "$fname" --quiet
    gunzip *.gz*
    ## Blast genome against PDV
    makeblastdb -in *.fna -dbtype nucl;
    blastn -task blastn -query PDV_HIMarthro.fasta -db *.fna  -outfmt '6 qaccver qlen saccver stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 8 -evalue 0.0001 -out outputblastn_${base_acc} ;
    ## delit gemone
    rm *.fna*;
    ## Filter ouputblast : evalue less than 0.0001 and lenght more than 165 pb
    #cat outputblastn_${base_acc} | awk -F "\t" '{if ($7>164 )print }'|awk -F "\t" '{if ($13<0.0001)print }' > outputblastn_filter_${base_acc}.txt;
    ## Add GenBank ID to blastn output
    find -iname "outputblastn_${base_acc}" -print0 | xargs -0 awk '
    FNR==1{sub(/\.txt$/, "", FILENAME)}
    FNR==1{sub("./", "", FILENAME)}
    FNR==1{sub("outputblastn_filter_", "", FILENAME)}
    {print $0 "\t" FILENAME}' > outputblastn_${base_acc}.txt;
    ## delete intermediate file
    rm outputblastn_${base_acc};
done

# Concatenate all outputblast
#cat outputblastn_filter_* > ALL_outputblastn_arthro_vs_PDV_${num}.txt
#rm outputblastn_filter_*
cat outputblastn_* > ALL_outputblastn_leftovers_${num}.txt
rm outputblastn_*