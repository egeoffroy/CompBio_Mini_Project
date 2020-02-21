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
    SRR1 = SRR[0:5] #first five values of SRR
    wget_command = 'wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/' + SRR1 + '/' + SRR+ '/'+ SRR+ '.sra'
    fastq_dump = 'fastq-dump -I --split-files ' + SRR + '.sra' #split the sra files into paired reads
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
    os.system(kallisto_build)
    kallisto = 'time kallisto quant -i HCMV_index.idx -o ./' + SRR+' -b 30 -t 4 '+ SRR + '_1.fastq ' + SRR+ '_2.fastq' #run kallisto for each of the SRRs
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
    original2 = open(SRR + '_2.fastq').readlines() #count the number of reads
    original = len(original1) + len(original2)
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
                output.write(str(file[i]) + '\n')
    output.close()

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
#            print(path)
            if len(path) == 0: #insert all kmer values as list and as set into stack
                for mer in mers:
                    substrings.append(([mer], set([mer]))) #add each mer to the list
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
def GenomeAssemblyShortest(sequences, final_sequence):
    if len(sequences) == 0:
        return final_sequence #return empty string
    elif (len(final_sequence) == 0):
        final_sequence = sequences.pop(0) #final sequences becomes the first sequence in the list of sequences                                                                                    
        return GenomeAssemblyShortest(sequences, final_sequence) #rerun with new final-sequence
    else:
        for i in range(len(sequences)):
            strings = sequences[i]
            string_length = len(strings) // 2 #length of strings divided by 2
            #print(string_length)
            for j in range(string_length): #for j in len of strings divided by 2
                c = len(strings) - j #c is length strings minus j
                if final_sequence.startswith(strings[j:]): #if the final sequence starts with a chunk of the next sequence
                    sequences.pop(i) #remove item from list
                    final_sequence = strings[:j] + final_sequence
                    return GenomeAssemblyShortest(sequences, final_sequence) #rerun function with string chunk added to final sequence
                #check from the opposite side
                if final_sequence.endswith(strings[:c]):
                    sequences.pop(i) #remove item from list
                    final_sequence = final_sequence + strings[c:]
                    return GenomeAssemblyShortest(sequences, final_sequence)



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
    outfile = open("assembly.fasta", 'a') #open the output file
    assembly = findAssembly(nodes_list) #find the Assembly for the Spades contigs
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
    outfile.close()

    with open('MiniProject.log', 'a') as output:
        output.write("There are " + str(count) + "  contigs > 1000 in the assembly" + '\n') #write out the number of contigs with bp >1000 to MiniProject.log
        output.close()
    return count

def contigs_length_count():
    file = open("./Spades/significant_contigs.txt").readlines()
    count=0
    for contig in file:
        number = len(contig)
        count+= number
    with open('MiniProject.log', 'a') as output:                                                                                                                                                  
        output.write("There are " + str(count) + " bp in the assembly" + '\n')  #count the number of bp in the assembly and write out to MiniProject.log
        output.close()

def blast():
    fasta = open("assembly.fasta").read()
    handle = NCBIWWW.qblast("blastn", "nr", fasta) #run blast against the assembled sequence
    with open("blast.xml", "w") as out_handle:
      out_handle.write(handle.read())
    out_handle.close()
    blast_qresult = SearchIO.read("blast.xml", "blast-xml")
    output = open('MiniProject', 'a')                                                                                                                                                         
    output.write('seq_title', 'align_len' , 'number_HSPs', 'topHSP_ident', 'topHSP_gaps', 'topHSP_bits', 'topHSP_expect')
    for i in range(0,10):
        hit = blast_qresult[i]
        blast_hsp = blast_qresult[i][0]
        output.write(str(hit.id) + ' ' + str(hit.seq_len) + ' '+ str(len(hit.hsps)) + ' ' + str(blast_hsp.ident_num) + ' ' + str(blast_hsp.gap_num) + ' ' + str(blast_hsp.bitscore)+ ' ' +str(blast_hsp.evalue))
                     #hit = blast_qresult[0]
    #blast_hsp = blast_qresult[0][0]
    #print(str(hit.id) + ' ' + str(hit.seq_len) + ' '+ str(len(hit.hsps)) + ' ' + str(blast_hsp.ident_num) + ' ' + str(blast_hsp.gap_num) + ' ' + str(blast_hsp.bitscore)+ ' '+str(blast_$    #for hit in blast_qresult[:10]: # id and sequence length of the first ten hits
    #  with open('MiniProject.log', 'a') as output:
    #      output.write("%s %i" %(hit.id, hit.seq_len, len(hit.hsps), hit.gap_num, hit.bitscore, hit.evalue) + '\n')
    output.close()






# ----------------- Theoretical Main Function ------------------                                                                                                                          
parser = argparse.ArgumentParser(description='Process some SRRs.')
parser.add_argument('SRRs', metavar='N', type=str, nargs='+', help='SRR values')
parser.add_argument('--download_files', help='Pull SRA Files from the database and split the paired reads')
args = parser.parse_args()
with open('MiniProject.log', 'a') as output:
  output.write('SRA values tested: ' + str(args.SRRs) + '\n')
  output.close()
#print(args.SRRs)

if args.download_files:
    for i in args.SRRs:
        Setup(i)
        with open('MiniProject.log', 'a') as output:
                output.write('SRA files downloaded')
                output.close()
                                                                                                                                                                                          
Transcriptome_index()
for i in args.SRRs:
    run_kallisto(i)

run_sleuth(args.SRRs)

for i in args.SRRs:
    bowtie_build(i)
    Count_bowtie(i)
    sam_to_fastq(i)

run_spades(args.SRRs[0], args.SRRs[1], args.SRRs[2], args.SRRs[3])                                                                                                                        
                     count = contigs_count()
contigs_length_count()
if count > 1: #if there is more than one significant contig, run the assembly
        run_find_Assembly()
else: #otherwise just make that one contig the assembly                                                                                                                                           
        file = open("./Spades/significant_contigs.txt").readlines()#
        for i in file:
                with open('assembly.fasta', 'a') as outfile:
                     outfile.write(i)
                     outfile.close()
#contigs_length_count()
blast()   
