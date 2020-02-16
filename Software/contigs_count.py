#Return number of contigs longer than 1000 bp
#Returns files of significant length contigs (>1000)
from Bio.Seq import Seq
from Bio import SeqIO
sequences = []
for record in SeqIO.parse("contigs.fasta", "fasta"):
    sequences.append(str(record.seq))

count = 0
outfile = open("./significant_contigs.txt", "w")
for contig in sequences:
    if len(contig) > 1000:
        count +=1
        outfile.write(contig + '\n')
print(count)
outfile.close()

