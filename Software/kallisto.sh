#!/bin/bash

time kallisto index -i HCMV_index.idx CDS_EF999921.txt

declare -a SRR
SRR=(SRR5660030 SRR5660033 SRR5660044 SRR5660045)
for i in "${SRR[@]}";
do
        time kallisto quant -i index.idx -o ./"${i}" -b 30 -t 4 ./ncbi_files/ncbi_files/"${i}"_1.fastq.gz ./ncbi_files/ncbi_files/"${i}"_2.fastq.gz
done

