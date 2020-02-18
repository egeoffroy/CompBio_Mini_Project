#!/bin/bash
#This file is deprecated and meant for downloading the entire SRA reads off of the NCBI SRA database
declare -a SRR
SRR=($1 $2 $3 $4)
for i in "${SRR[@]}"
do
  v2=${i::-4}
  wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${v2}/${i}/${i}.sra

done
#wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR566/SRR5660030/SRR5660030.sra
#wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR566/SRR5660033/SRR5660033.sra
#wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR566/SRR5660044/SRR5660044.sra
#wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR566/SRR5660045/SRR5660045.sra
#curl -O http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRR5660030&format=fastq
#curl -O http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRR5660033&format=fastq
#curl -O http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRR5660044&format=fastq
#curl -O http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRR5660045&format=fastq

curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=EF999921&rettype=gb&retmode=txt">E$
