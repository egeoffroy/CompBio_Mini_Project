#!/bin/bash

#This is an extra step in case the user only has the SRA files and needs to separate them into paired reads
declare -a SRRs
SRRs=(SRR5660030 SRR5660033 SRR5660044 SRR5660045)
for i in "${SRRs[@]}":
do
        fastq-dump --split-files -O ./ "${i}"
        echo "$i";
done
