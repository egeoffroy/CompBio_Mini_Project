#!/bin/bash
#cd ..
declare -a SRR
SRR=($1 $2 $3 $4)
bowtie2-build ./EF999921.fasta EF999921
for i in "${SRR[@]}"
do
        bowtie2 --quiet -x EF999921 -1 ./${i}_1.fastq.gz -2 ./${i}_2.fastq.gz -S EF999921_${i}.sam
done
