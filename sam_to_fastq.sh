#!/bin/bash

#This file takes the sam output files from Bowtie2 and converts them into fastq files for Spades

sam=$1
cat EF999921_${1}.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > EF999921_${1}.fastq
