#Create Transcriptome Index for HCMV
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "ecgeoffroy@gmail.com"

handle = Entrez.efetch(db="nucleotide", id=["EF999921"], rettype="fasta") #obtain plain text records of the GenBank ids in FASTA format from NCBI's [Nucleotide] database.
records = list(SeqIO.parse(handle, "fasta"))
SeqIO.write(records[0],'EF999921.fasta', 'fasta') #write out the fasta file for bowtie2

outfile = open("CDS_EF999921.txt", 'w') #open write out file
i = 0 #use i to count the number of CDS entries written out to the file
for rec in SeqIO.parse('EF999921.gb', 'genbank'): #open the genBank file
        if rec.features: #look at the genbank features
#Create Transcriptome Index for HCMV
from Bio import SeqIO
from Bio import Entrez
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
