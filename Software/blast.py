from Bio.Blast import NCBIWWW
fasta = open("assembly.fasta").read()
handle = NCBIWWW.qblast("blastn", "nr", fasta)

with open("blast.xml", "w") as out_handle:
  out_handle.write(handle.read())
out_handle.close()

from Bio import SearchIO
blast_qresult = SearchIO.read("blast.xml", "blast-xml")
for hit in blast_qresult[:10]:   # id and sequence length of the first five hits
  print("%s %i" % (hit.id, hit.seq_len, len(hit.hsps), hit.ident_num, hit.gap_num, hit.bitscore, hit.evalue))
