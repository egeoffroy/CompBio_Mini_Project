#!/bin/bash

declare -a SRRs
SRRs=($1 $2 $3 $4)
echo $SRRs
cd ..
log=miniProject.log
#curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=EF999921&rettype=gb&retmode=txt">EF999921.gb
#bash ./Software/download_files.sh #deprecated line
var=$(python ./Software/Transcriptome_index.py 2>&1)
echo 'The HCMV genome (EF999921) has ' $var >> $log
#bash ./Software/split_paired.sh
gzip SRR*
bash ./Software/kallisto.sh $1 $2 $3 $4
bash ./Software/bowtie2.sh $1 $2 $3 $4

#For each of the SRR values, determine which donor it is from and find the number of read pairs before and after Bowtie2 filtering
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


#spades -k 55,77,99,127 -t 2 --pe1-1 ./ncbi_files/ncbi_files/SRR5660030_1.fastq.gz --pe1-2 ./ncbi_files/ncbi_files/SRR5660030_2.fastq.gz  --pe2-1 ./ncbi_files/ncbi_files/SRR5660033_1.fastq.gz --pe2-2 ./ncbi_files/ncbi_files/SRR5660033_2.fastq.gz --pe3-1 ./ncbi_files/ncbi_files/SRR5660044_1.fastq.gz --pe3-2 ./ncbi_files/ncbi_files/SRR5660044_2.fastq.gz --pe4-1 ./ncbi_files/ncbi_files/SRR5660045_1.fastq.gz --pe4-2 ./ncbi_files/ncbi_files/SRR5660045_2.fastq.gz -o ./
#spades -k 55,77,99,127 -t 2 --pe1-1 SRR5660030_1.fastq.gz -pe1-2 SRR5660030_2.fastq.gz --pe2-1 SRR5660033_1.fastq.gz -pe2-2 SRR5660033_2.fastq.gz --pe3-1 SRR5660044_1.fastq.gz --pe3-2 SRR5660044_2.fastq.gz --pe4-1 SRR5660045_1.fastq.gz --pe4-2 SRR5660045_2.fastq.gz -o ./
echo 'There are __  contigs > 1000 bp in the assembly ' >> $log
echo 'There are __  bp in the assembly' >> $log
