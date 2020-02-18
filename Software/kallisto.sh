#!/bin/bash

#Create the transcriptome reference index
time kallisto index -i HCMV_index.idx EF999921.fasta

# For each of the user inputted SRRs, run kallisto quantification with their read paired separated files
declare -a SRR
SRR=($1 $2 $3 $4)
for i in "${SRR[@]}";
do
        time kallisto quant -i HCMV_index.idx -o ./"${i}" -b 30 -t 4 ./"${i}"_1.fastq.gz ./"${i}"_2.fastq.gz
done
