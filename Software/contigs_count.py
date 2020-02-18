#Return number of contigs longer than 1000 bp
#Returns files of significant length contigs (>1000)

#Import libraries
from Bio.Seq import Seq
from Bio import SeqIO

sequences = [] #make empty array
for record in SeqIO.parse("contigs.fasta", "fasta"):
    sequences.append(str(record.seq)) #add fasta sequences to the array sequences

count = 0
outfile = open("./significant_contigs.txt", "w") # open new file
for contig in sequences: #for each of the sequences (contigs)
    if len(contig) > 1000: #check if the length of bp is > 1000
        count +=1 #if it is, add one to count and write out the contig to the file
        outfile.write(contig + '\n')
print(count) #print out the count of contigs > 1000bp
outfile.close()

