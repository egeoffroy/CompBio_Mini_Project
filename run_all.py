import os
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq

path = os.getcwd()

def Setup(SRR):
    #SRR1 = SRR[0:5] #first five values of SRR --> deprecated
    #https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1
    wget_command = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR+ '/'+ SRR+ '.1 '
    fastq_dump = 'fastq-dump -I --split-files ' + SRR + '.1' #split the sra files into paired reads
    os.system(wget_command) #run the commands
    os.system(fastq_dump)

def Transcriptome_index():
    fasta = 'curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=EF999921&rettype=gb&retmode=txt">EF999921.gb'
    os.system(fasta) #pull the genbank record for EF999921
    outfile = open("CDS_EF999921.fasta", 'w') #open the output file
    i = 0
    Entrez.email = "ecgeoffroy@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=["EF999921"], rettype="fasta") #obtain plain text records of the GenBank ids in FASTA format from NCBI's [Nucleotide] database.
    records = list(SeqIO.parse(handle, "fasta"))
    fasta_out = open("EF999921.fasta", 'w')
    fasta_out.write(str(records[0].description) + '\n' + str(records[0].seq)) #write out the fasta sequence to the file
    fasta_out.close()
    SeqIO.write(records[0],'EF999921.fasta', 'fasta')
    for rec in SeqIO.parse("EF999921.gb", "genbank"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS": #find the CDS regions and add them to the transcriptome index file
                    outfile.write('>' + str(feature.qualifiers['protein_id']).replace('[', '').replace(']', '').replace("'", '') + '\n' + str(feature.location.extract(rec).seq) + '\n') 
                    i+=1

    with open('MiniProject.log', 'a') as output:
        output.write('The HCMV genome (EF999921) has '+ str(i) + '\n') #output how many CDS regions there are in the EF999921 genbank file
        output.close()
    outfile.close()

def run_kallisto(SRR):
    kallisto_build = 'time kallisto index -i HCMV_index.idx CDS_EF999921.fasta' ##create the transcriptome index in kallisto
    os.system(kallisto_build) #run command 
    kallisto = 'time kallisto quant -i HCMV_index.idx -o ./' + SRR+' -b 30 -t 4 '+ SRR + '.1_1.fastq ' + SRR+ '.1_2.fastq' #run kallisto for each of the SRRs
    os.system(kallisto)

def sam_to_fastq(SRR):
    fastq = 'bash sam_to_fastq.sh ' + SRR #convert the .sam bowtie output files to .fastq
    os.system(fastq)

def run_spades(SRR1, SRR2, SRR3, SRR4):
    spades_command = 'spades -k 33,55,77,99,127 -t 2 -s EF999921_'+ SRR1 + '.fastq -s EF999921_'+ SRR2 + '.fastq -s EF999921_' + SRR3 + '.fastq -s EF999921_' + SRR4 + '.fastq -o Spades/'    
    with open('MiniProject.log','a') as output: #run SPades and print out the command to MiniProject.log
        output.write(spades_command + '\n')
        output.close()
    os.system(spades_command)

def bowtie_build(SRR):
    #builds initial index using reference genome-- use transcriptome_index CDS fasta file?
    bowtie_command = 'bowtie2-build ./CDS_EF999921.fasta EF999921'
    os.system(bowtie_command)
    #maps transcriptome reads to the index we just created, generating sam file
    bowtie_command2 = 'bowtie2 --quiet --no-unal -x EF999921 -1 '+ SRR+ '.1_1.fastq -2 ' + SRR+ '.1_2.fastq -S EF999921_' + SRR+ '.sam'
    os.system(bowtie_command2)

def Count_bowtie(SRR):
    bowtie_SRR = open('EF999921_' + SRR + '.fastq').readlines()
    donor = ''
    if SRR == 'SRR5660030':
        donor += 'Donor 1 (2dpi)' #This is specific for the particular reads --> probably could make each run one of the values 
    elif SRR == 'SRR5660033':
        donor += 'Donor 1 (6dpi)'
    elif SRR == 'SRR5660044':
        donor += 'Donor 3 (2dpi)'
    elif SRR == 'SRR5660045':
        donor += 'Donor 3 (6dpi)'
    len_bowtie = (len(bowtie_SRR)/4)
    original1 = open(SRR + '.1_1.fastq').readlines()                                                                                                                                            
    original2 = open(SRR + '.1_2.fastq').readlines() #count the number of reads
    original = (len(original1) + len(original2))/4
    #write out to the log file
    with open('MiniProject.log', 'a') as output:
        output.write(donor + " had " + str(original) + ' read pairs before Bowtie2 filtering and ' + str(len_bowtie) + ' read pairs after \n')
        output.close()
        
def run_sleuth(SRR):
    run_R = 'Rscript make_sample_covariates.R ' + SRR[0] + ' ' + SRR[1]+ ' ' + SRR[2] + ' ' + SRR[3] #make the sample_covariance.txt file
    os.system(run_R)
    run_more_R = 'Rscript sleuth.R' #run sleuth
    os.system(run_more_R)                                                                                                                                                                     
    file = open('sleuth_output.txt').readlines()
    with open('MiniProject.log' ,'a') as output:
        for i in range(len(file)):
                output.write(str(file[i]) + '\n') #print out the top sleuth hits
    output.close()

#Assemble contigs > 1000 bp
def concat_contigs():
    file = open('./Spades/significant_contigs.txt') #input file                                                             
    file = file.read().split() #split the file
    final_seq = '' #create empty string for the assembly to be inputted to
    for i in range(len(file)):
        final_seq += file[i] + 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    output1 = open('assembly.fasta' , 'a')                                #write out the final sequence to the assembly.fasta file                                                  
    output1.write(final_seq)
    output1.close()
    
def contigs_count():
    sequences = [] #make empty array
    for record in SeqIO.parse("./Spades/contigs.fasta", "fasta"):
        sequences.append(str(record.seq)) #add fasta sequences to the array sequences
    count = 0
    outfile = open("./Spades/significant_contigs.txt", "w") # open new file
    for contig in sequences: #for each of the sequences (contigs)
        if len(contig) > 1000: #check if the length of bp is > 1000
            count +=1 #if it is, add one to count and write out the contig to the file                                                                                                                
            outfile.write(contig + '\n')
    outfile.close()

    with open('MiniProject.log', 'a') as output:
        output.write("There are " + str(count) + "  contigs > 1000 in the assembly" + '\n') #write out the number of contigs with bp >1000 to MiniProject.log
        output.close()
    return count

def contigs_length_count():
    file = open("./Spades/significant_contigs.txt").readlines()
    count=0
    for contig in file:
        number = len(contig) #add up the number of bp in the assembly
        count+= number
    with open('MiniProject.log', 'a') as output:                                                                                                                                                  
        output.write("There are " + str(count) + " bp in the assembly" + '\n')  #count the number of bp in the assembly and write out to MiniProject.log
        output.close()

def blast():
    fasta = open("assembly.fasta").read()
    handle = NCBIWWW.qblast("blastn", "nr", fasta, entrez_query='10292[taxid]') #run blast against the assembled sequence
    with open("blast.xml", "w") as out_handle:
      out_handle.write(handle.read())
    out_handle.close()
    blast_qresult = SearchIO.read("blast.xml", "blast-xml")
    output = open('MiniProject.log', 'a')                                                                                                                                                         
    output.write('seq_title'+ 'align_len' + 'number_HSPs'+ 'topHSP_ident'+ 'topHSP_gaps'+ 'topHSP_bits'+ 'topHSP_expect')
    max_blast_id = 10
    if len(blast_qresult) < 10: #prevents program from crashing when there are less than 10 results
        max_blast_id = len(blast_qresult)
    for i in range(0,max_blast_id):
        hit = blast_qresult[i]
        blast_hsp = blast_qresult[i][0]
        output.write(str(hit.id) + ' ' + str(hit.seq_len) + ' '+ str(len(hit.hsps)) + ' ' + str(blast_hsp.ident_num) + ' ' + str(blast_hsp.gap_num) + ' ' + str(blast_hsp.bitscore)+ ' ' +str(blast_hsp.evalue))
    output.close()






# ----------------- Theoretical Main Function ------------------                                                                                                                          
parser = argparse.ArgumentParser(description='Process some SRRs.')
parser.add_argument('SRRs', metavar='N', type=str, nargs='+', help='SRR values')
parser.add_argument('--download_files', default='N', help='Pull SRA Files from the database and split the paired reads')
args = parser.parse_args()

# Call all of the functions using the user given SRRs
with open('MiniProject.log', 'a') as output:
  output.write('SRA values tested: ' + str(args.SRRs) + '\n')
  output.close()

if args.download_files != 'N':
    for i in args.SRRs:
        Setup(i)
        with open('MiniProject.log', 'a') as output:
                output.write(i + ' SRA file downloaded \n')
                output.close()
                                                                                                                                                                                          
Transcriptome_index()
for i in args.SRRs:
    run_kallisto(i)

run_sleuth(args.SRRs)
for i in args.SRRs:
    bowtie_build(i)
    sam_to_fastq(i)
    Count_bowtie(i)
    

run_spades(args.SRRs[0], args.SRRs[1], args.SRRs[2], args.SRRs[3])                                                                                                                        
count = contigs_count()
contigs_length_count()
if count > 1: #if there is more than one significant contig, run the assembly
        concat_contigs()
else: #otherwise just make that one contig the assembly                                                                                                                                           
        file = open("./Spades/significant_contigs.txt").readlines()#
        for i in file:
                with open('assembly.fasta', 'a') as outfile:
                     outfile.write(i)
                     outfile.close()
blast()   
