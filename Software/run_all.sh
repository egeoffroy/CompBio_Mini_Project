#!/bin/bash

declare -a SRRs
SRRs=($1 $2 $3 $4)
echo $SRRs
cd ..
log=miniProject.log
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=EF999921&rettype=gb&retmode=txt">EF999921.gb
var=$(python ./Software/Transcriptome_index.py 2>&1)
echo 'The HCMV genome (EF999921) has ' $var >> $log
#bash ./Software/split_paired.sh
#gzip SRR*
bash ./Software/kallisto.sh $1 $2 $3 $4
echo $PWD
Rscript ./Software/make_sample_covariates.R $1 $2 $3 $4
Rscript ./Software/sleuth.R
bash ./Software/bowtie2.sh $1 $2 $3 $4

#For each of the SRR values, determine which donor it is from and find the number of read pairs before and after Bowtie2
for i in "${SRRs[@]}";
do
        if [ "${i}" == 'SRR5660030' ]
        then
                donor=('Donor 1 (2dpi)')
        elif [ "${i}" == 'SRR5660033' ]
        then
                donor=('Donor 1 (6dpi)')
        elif [ "${i}" == 'SRR5660044' ]
        then
                donor=('Donor 3 (2dpi)')
        elif [ "${i}" == 'SRR5660045' ]
        then
                donor=('Donor 3 (6dpi)')
        fi
        after1=$(wc -l < EF999921_${i}.sam)
        after=$((after1 - 3))
        before1=$(wc -l < ./"${i}"_1.fastq.gz)
        before2=$(wc -l < ./"${i}"_2.fastq.gz)
        #before1=$((before1 / 2))
        #before2=$((before2 / 2))
        before=$((before1 + before2))

        echo ${donor} 'had' ${before} 'read pairs before Bowtie2 filtering and' ${after} 'read pairs after' >> $log
done

cat EF99921_SRR5660030.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > EF99921_SRR5660030.fastq
cat EF99921_SRR5660033.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > EF99921_SRR5660033.fastq
cat EF99921_SRR5660044.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > EF99921_SRR5660044.fastq
cat EF99921_SRR5660045.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > EF99921_SRR5660045.fastq


spades -k 55,77,99,127 -t 2 -1 EF99921_SRR5660030.fastq -2 EF999921_SRR5660033.fastq -3 EF999921_SRR5660044.fastq -4 EF999921_SRR5660045.fastq -o ./
echo 'spades -k 55,77,99,127 -t 2 -1 EF99921_SRR5660030.fastq -2 EF999921_SRR5660033.fastq -3 EF999921_SRR5660044.fastq -4 EF999921_SRR5660045.fastq -o ./' >> $log

var1=$(python ./Software/contigs_count.py 2>&1)
echo 'There are' ${var1}'contigs > 1000 bp in the assembly ' >> $log
var2=$(python ./Software/contigs_length_count.py 2>&1)
echo 'There are' ${var2} 'bp in the assembly' >> $log

blast=$python( ./Software/blast.py 2>&1) >> $log
