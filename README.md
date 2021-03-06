# Computional Biology Mini Project

## Installation
To run this program from your own working directory, clone this repository to your workspace. 

```
  git clone https://github.com/egeoffroy/CompBio_Mini_Project.git
```

### Software:
#### Python
##### os -- https://docs.python.org/3/library/os.html
##### argparse -- https://docs.python.org/3/library/argparse.html
##### Biopython -- http://biopython.org/DIST/docs/tutorial/Tutorial.html
###### The packages within Biopython that are used within the program are SeqIO, SearchIO, NCBIWWW, Seq, and Entrez. 
#### Kallisto -- https://pachterlab.github.io/kallisto/
#### Bowtie2 -- http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
#### SPades -- http://cab.spbu.ru/software/spades/

## Directions
To run this program, move the SRR files you would like to run the program with in the directory CompBio_Mini_Project. This can be done by running:
```
  mv SRR_file CompBio_Mini_Project/SRR_file
```
From there, set your working directory to CompBio_Mini_Project and call the function. The program requires four SRR values as input. 
```
  cd CompBio_Mini_Project
  python run_all.py SRR1 SRR2 SRR3 SRR4
```

Or alternatively, you can use the flag --download_files to pull the SRA files from the database. 
```
  cd CompBio_Mini_Project
  python run_all.py SRR1 SRR2 SRR3 SRR4 --download_files Y
```

Or to run it from the background,
```
  cd CompBio_Mini_Project
  nohup python run_all.py SRR1 SRR2 SRR3 SRR4 --download_files Y &
```
## Test Data 
To run the software with the test data provided:
```
  cd CompBio_Mini_Project
  python run_all.py SRR5660030 SRR5660033 SRR5660044 SRR5660045
```

This test data consists of the first 50000 lines for the four SRR paired read files. This means there are 12500 read pairs for each SRR. When running with this tool, expect only one significant contig with length > 1000 bp after running SPades. 

## Scripts

#### run_all.py : a python script that calls the scripts make_sample_covariates.R, sleuth.R, and sam_to_fastq.sh. It also includes functions that download various SRA files, pulls fasta files and GenBank entries, runs kallisto, runs Bowtie2, run SPades, assembles the reads, counts the number of reads and number of base pairs in the assembly, and runs blast with the assembly against the nr database. 

#### make_sample_covariates.R : a R script that makes a sleuth prep file with columns sample, path, and condition. Sample is the name of the SRR. Path is the file path to the directory of the kallisto output, and condition identifies whether the read is from day 2 or day 6

#### sleuth.R : a R script that creates a sleuth object and identifies differentially expressed genes between two timepoints.




## Output

#### MiniProject.log : a log file that contains various information from the tests with your inputted data. This includes the number of base pairs in the assembly, number of contigs with bp > 1000, top ten sleuth hits, top ten BLAST hits, etc. 

 EF999921.gb : a genbank file for EF999921
 
 CDS_EF999921.fasta : a fasta file for only the coding regions of EF999921 
 
 sample_covariates.txt : a file containing the sample names, paths, and conditions in order to run sleuth.R
 
 sleuth_output.txt : a file containing the significant hits after running sleuth. 

 significant_contigs.txt : a file listing the contigs with more than 1000 bp. 
 
 assembly.fasta : the assembled significant contigs. 

 blast.xlm : the blast output file.










### 15pts
### Due via Sakai by 11:59PM 2/27
### Description: The mini-project will focus on integrating some of the software tools discussed in class. Students will work
independently to develop a Python wrapper to automate the execution of the software tools. You will create a GitHub
repo and post your code there (Introduction to GitHub presented on 2/13 in class). You should include straight-forward
documentation and sample data in your repository. The code should have adequate documentation in the repository’s
README.md file such that anyone could download the code, install it or get it to run, and run through the test data
without encountering any problems; this should be thought of as a User’s Manual. If I cannot run the code with the test
data, I cannot grade the functionality of the code. This means that you cannot hardcode paths in your code!
The code will be graded as follows:
• 3 pts. Code comments
• 3 pts. Documentation -- including what tools need to be installed and how to use the code
• 1 pts. Test data (include a small subset of input reads so I can test your code quickly, but include enough reads so the
entire pipeline runs. Note, files in your GitHub repo must be less than 50MB.)
• 8 pts. Functionality (include miniProject.log in your GitHub repo with the requested output from running your pipeline
with all input reads, so I can check your answers.)
Submit the URL to your GitHub repo to Sakai.
Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV.
From Wikipedia: Although they may be found throughout the body, HCMV infections are frequently associated with the
salivary glands. HCMV infection is typically unnoticed in healthy people, but can be life-threatening for the
immunocompromised, such as HIV-infected persons, organ transplant recipients, or newborn infants. Congenital
cytomegalovirus infection can lead to significant morbidity and even death. After infection, HCMV remains latent within
the body throughout life and can be reactivated at any time. Eventually, it may cause mucoepidermoid carcinoma and
possibly other malignancies such as prostate cancer.
Cheng et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/29158406) produced the transcriptomes of HCMV post
infection. Write a Python script to automate the following and produce the output file requested named
“miniProject.log” and other output files in a folder named “miniProject_FirstName_LastName” [where you’ve indicated
your first and last name]. ALL results generated by you or programs called should be written to this folder. The easiest
way to guarantee this is to create the folder and then move into it via an os.system call using cd (change directory).

# Questions
1. We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi). First, retrieve the following
transcriptomes from two patient donors from SRA and convert to paired-end fastq files. You can use wget (by
constructing the path based on the SRR numbers for each of these samples).
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363
Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374
Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375
2. We will quantify TPM in each sample using kallisto, but first we need to build a transcriptome index for HCMV (NCBI
accession EF999921). Use Biopython to retrieve and generate the appropriate input and then build the index with
kallisto (https://pachterlab.github.io/kallisto/). You will need to extract the CDS features from the GenBank format.
Write the following to your log file (replace # with the number of coding sequences in the HCMV genome):
2
The HCMV genome (EF99921) has # CDS.
3. Quantify the TPM of each CDS in each transcriptome using kallisto and use these results as input to find differentially
expressed genes between the two timepoints (2pi and 6dpi) using the R package sleuth
(https://pachterlab.github.io/sleuth/about). Write the following details for each significant transcript (FDR < 0.05) to
your log file, include a header row, and tab-delimit each item:
target_id test_stat pval qval
4. It has been proposed that HCMV disease and pathogenesis may be related to the genetic diversity of the virus
(Renzette et al. https://www.ncbi.nlm.nih.gov/pubmed/25154343/). Which publicly available strains are most similar to
these patient samples? To compare to other strains, we will assemble these transcriptome reads. We don’t expect
assembly to produce the entire genome, but enough to be useful in BLAST. Virus sequencing experiments often include
host DNAs. It is difficult to isolate the RNA of just the virus (as it only transcribes during infection of the host cell). Before
assembly, let’s make sure our reads map to HCMV. Using Bowtie2, create an index for HCMV (NCBI accession EF999921).
Next, save the reads that map to the HCMV index for use in assembly. Write to your log file the number of reads in each
transcriptome before and after the Bowtie2 mapping. For instance, if I was looking at the Donor 1 (2dpi) sample, I would
write to the log (numbers here are arbitrary):
Donor 1 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read
pairs after.
5. Using the Bowtie2 output reads, assemble all four transcriptomes together to produce 1 assembly via SPAdes. Write
the SPAdes command used to the log file.
6. Write Python code to calculate the number of contigs with a length > 1000 and write the # out to the log file as
follows:
There are # contigs > 1000 bp in the assembly.
For steps 7-9, you will only consider those contigs > 1000 bp in length.
7. Write Python code to calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in
length) and write this # out to the log file as follows:
There are # bp in the assembly.
8. Write Python code to concatenate all of the contigs > 1000 bp in length into 1 fasta sequence. Separate each contig by
a stretch of 50 N’s.
9. Using this concatenated fasta file, blast via NCBIWWW.qblast to query the nr database limited to members of the
Herpesviridae family to identify the top 10 hits. For the top 10 hits, write the following to your log file: sequence title,
alignment length, number of HSPs, and for the top HSP: HSP identities, HSP gaps, HSP bits, and HSP expect scores.
Include the following header row and tab-delimit each item:
seq_title align_len number_HSPs topHSP_ident topHSP_gaps topHSP_bits topHSP_expect
