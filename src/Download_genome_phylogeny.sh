acc_file=list_species_phylo.txt

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
    rm *.gz*
done
