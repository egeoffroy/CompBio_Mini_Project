import os
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
#import csv
path = os.getcwd()

def Setup(SRR):
    SRR1 = SRR[0:5]
    wget_command = 'wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/' + SRR1 + '/' + SRR+ '/'+ SRR+ '.sra'
    fastq_dump = 'fastq-dump -I --split-files ' + SRR + '.sra'
    os.system(wget_command)
    os.system(fastq_dump)

def Transcriptome_index():
    fasta = 'curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=EF999921&rettype=gb&retmode=txt">EF999921.gb'
    os.system(fasta)

    Entrez.email = "ecgeoffroy@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=["EF999921"], rettype="fasta") #obtain plain text records of the GenBank ids in FASTA format from NCBI's [Nucleotide] database.
    records = list(SeqIO.parse(handle, "fasta"))
    SeqIO.write(records[0],'EF999921.fasta', 'fasta') #write out the fasta file for bowtie2
    outfile = open("CDS_EF999921.txt", 'w') #open write out file
    i = 0 #use i to count the number of CDS entries written out to the file
    for rec in SeqIO.parse('EF999921.gb', 'genbank'): #open the genBank file
            if rec.features: #look at the genbank features
                    for feature in rec.features:
                            if feature.type == "CDS": #only write out the CDS features
                                    outfile.write(str(feature))
                                    i+=1
    print(i)
    outfile.close()


def run_kallisto(SRR):
    kallisto_build = 'time kallisto index -i HCMV_index.idx EF999921.fasta'
    os.system(kallisto_build)
    kallisto = 'time kallisto quant -i HCMV_index.idx -o ./' + SRR+' -b 30 -t 4 '+ SRR + '_1.fastq ' + SRR+ '_2.fastq'
    os.system(kallisto)

def sam_to_fastq(SRR):
    fastq = 'bash sam_to_fastq.sh ' + SRR
    os.system(fastq)

def run_spades(SRR1, SRR2, SRR3, SRR4):
    spades_command = 'spades -k 33,55,77,99,127 -t 2 -s EF999921_'+ SRR1 + '.fastq -s EF999921_'+ SRR2 + '.fastq -s EF999921_' + SRR3 + '.fastq -s EF999921_' + SRR4 + '.fastq -o Spades/'
    with open('MiniProject.log','a') as output:
        output.write(spades_command)
        output.close()
    os.system(spades_command)
    #writes spades command to log file and executes

def bowtie_build(SRR):
    #builds initial index using reference genome
    bowtie_command = 'bowtie2-build ./EF999921.fasta EF999921'
    os.system(bowtie_command)
    #maps transcriptome reads to the index we just created, generating sam file
    os.system('fastq-dump -I --split-files' + SRR + '.sra')
    bowtie_command2 = 'bowtie2 --quiet -x EF999921 -1 '+ SRR+ '_1.fastq -2 ' + SRR+ '_2.fastq -S EF999921_' + SRR+ '.sam'
    os.system(bowtie_command2)

def Count_bowtie(SRR):
    bowtie_SRR = open('EF999921_' + SRR + '.sam').readlines()
    donor = ''
    if SRR == 'SRR5660030':
        donor += 'Donor 1 (2dpi)'
    elif SRR == 'SRR5660033':
        donor += 'Donor 1 (6dpi)'
    elif SRR == 'SRR5660044':
        donor += 'Donor 3 (2dpi)'
    elif SRR == 'SRR5660045':
        donor += 'Donor 3 (6dpi)'
    len_bowtie = len(bowtie_SRR)
    original1 = open(SRR + '_1.fastq').readlines()
    original2 = open(SRR + '_2.fastq').readlines()
    original = len(original1) + len(original2)
    #write out to the log file
    with open('MiniProject.log', 'a') as output:
        output.write(donor + "had" + str(original) + 'read pairs before Bowtie2 filtering and ' + str(len_bowtie) + 'read pairs after \n')
        output.close()

def run_sleuth(SRR):
    run_R = 'Rscript make_sample_covariates.R ' + SRR[0] + ' ' + SRR[1]+ ' ' + SRR[2] + ' ' + SRR[3]
    os.system(run_R)
    run_more_R = 'Rscript sleuth.R'
    os.system(run_more_R)
 
def findAssembly(nodes_list):
    for k in reversed(range(2, len(nodes_list[0])+1)): #from k starting at len of list +1 to 2 bc of the function reverse (account for index start at 0)
        mers = set()#create empty set
        for s in nodes_list: #for value in nodes_list
           for i in range(len(s)-k+1): #for i in length of value in list minus k+1
                mers.add(s[i:i+k]) #add kmer to set

        final_string = None
        substrings = []
        substrings.append(([], set()))

        while len(substrings) > 0 and final_string == None:
            path, vs = substrings.pop() #assign values
            print(path)
            if len(path) == 0: #insert all kmer values as list and as set into stack
                for mer in mers:
                    substrings.append(([mer], set([mer])))  #add each mer to the list
            else:
                mer = path[-1] #mer is last value in list path
                print(mer)
                for a in "ACGT": #for each nucleotide
                    nmer = mer[1:] + a #next mer is the mer[1:] and the next nucleotide a
                    if nmer in mers and nmer != str(Seq(mer).reverse_complement()): #if next mer is in mers and it is not equal to a reverse complement
                        if nmer == path[0]: #if next mer is equal to
                            final_string = list(path)
                            break
                        elif not nmer in vs:
                            substrings.append((path + [nmer], vs.union(set([nmer]))))
        #print(stack)
        if final_string != None: #if final string isn't empty
            output = final_string[0] #assign first value to output
            for i in range(1, len(final_string)):
                output += final_string[i][-1]  #add each value in final_string to the output but subtract the first letter
            #print(output)
            output = output[:-(k-1)] #remove
            doutput = output + output
            result = True #create boolean variable

            for dna in nodes_list:
                if doutput.find(dna) == -1 and doutput.find(str(Seq(dna).reverse_complement())) == -1:
                    result = False
            if result:
                return(output)
                break
 def run_find_Assembly():
    file = open("./Spades/significant_contigs.txt") #input file
    file = file.read().split()

    src = [] #set of reverse complements
    for i in file:
        sequence = Seq(i)
        reverse = sequence.reverse_complement()
        src.append(str(reverse))

    a = set(file)
    nodes = set(src).union(set(a))
    nodes_list = list(nodes)
    outfile = open("assembly.fasta", 'w') #open the output file
    assembly = findAssembly(nodes_list)
    for sub in range(0, len(assembly), 50):
        outfile.write(assembly[sub:sub+50]) #write out the file
    outfile.close()

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
    with open('MiniProject.log', 'a') as output:
        output.write("There are %d contigs > 1000 in the assembly"%(count) + '\n')
    #print(count) #print out the count of contigs > 1000bp
    outfile.close()

def contigs_length_count():
    file = open("./Spades/significant_contigs.txt").readlines()
    count=0
    for contig in file:
        number = len(contig)
        count+= number
    with open('MiniProject.log', 'a') as output:
        output.write("There are %d bp in the assembly"%(count) + '\n')

def blast():
    fasta = open("assembly.fasta").read()
    handle = NCBIWWW.qblast("blastn", "nr", fasta)
    with open("blast.xml", "w") as out_handle:
      out_handle.write(handle.read())
    out_handle.close()
    blast_qresult = SearchIO.read("blast.xml", "blast-xml")
    for hit in blast_qresult[:10]:   # id and sequence length of the first five hits
      with open('MiniProject.log', 'a') as output:
          output.write("%s %i" %(hit.id[0:10], hit.seq_len[0:10], len(hit.hsps)[0:10], hit.ident_num[0:10], hit.gap_num[0:10], hit.bitscore[0:10], hit.evalue[0:10]) + '\n')
          output.close()








parser = argparse.ArgumentParser(description='Process some SRRs.')
parser.add_argument('SRRs', metavar='N', type=str, nargs='+', help='SRR values')
parser.add_argument('--download_files', help='Pull SRA Files from the database and split the paired reads')

args = parser.parse_args()
with open('MiniProject.log', 'a') as output:
  output.write('SRA values tested: ' + str(args.SRRs) + '\n')
print(args.SRRs)
#SRR = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

if args.download_files:
    for i in args.SRRs:
        Setup(i)

Transcriptome_index()
for i in args.SRRs:
    run_kallisto(i)

run_sleuth(SRR)
for i in args.SRRs:
    bowtie_build(i)
    Count_bowtie(i)
    sam_to_fastq(i)

run_spades(args.SRRs[0], args.SRRs[1], args.SRRs[2], args.SRRs[3])
contigs_count()
run_find_Assembly()
blast()
