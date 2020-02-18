#!/bin/bash

#This is an extra step in case the user only has the SRA files and needs to separate them into paired reads
declare -a SRRs
SRRs=($1 $2 $3 $4)
for i in "${SRRs[@]}":
do
        fastq-dump --split-files -O ./ "${i}"
        echo "$i";
done
