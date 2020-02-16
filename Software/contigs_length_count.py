# Count length of each contig > 1000 bp

file = open("significant_contigs.txt").readlines()
count=0
for contig in file:
    number = len(contig)
    count+= number
print(count)
