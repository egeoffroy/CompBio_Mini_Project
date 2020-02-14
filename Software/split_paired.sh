#!/bin/bash
declare -a SRRs
SRRs=(SRR5660030 SRR5660033 SRR5660044 SRR5660045)
for i in "${SRRs[@]}":
do
        fastq-dump --split-files -O ./ "${i}"
        echo "$i";
done
