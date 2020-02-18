#!/bin/bash
#A script to run Bowtie2 with the various reads 
declare -a SRR
SRR=($1 $2 $3 $4)
bowtie2-build ./EF999921.fasta EF999921 #run build 
#for each of the SRRs, run Bowtie --> using single reads instead of paired reads
for i in "${SRR[@]}" 
do
        bowtie2 --quiet -x EF999921 -1 ./${i}_1.fastq -2 ./${i}_2.fastq -S EF999921_${i}.sam #output is a .sam file
done
