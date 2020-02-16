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
            print(output) 
            output = output[:-(k-1)] #remove 
            doutput = output + output
            result = True #create boolean variable
            
            for dna in nodes_list:
                if doutput.find(dna) == -1 and doutput.find(str(Seq(dna).reverse_complement())) == -1:
                    result = False
            if result:
                return(output)
                break
    

file = open("./significant_contigs.txt") #input file
file = file.read().split()
        
src = [] #set of reverse complements
for i in file:
    sequence = Seq(i)
    reverse = sequence.reverse_complement()
    src.append(str(reverse))
    
a = set(file)
nodes = set(src).union(set(a))
nodes_list = list(nodes)
outfile = open("assembly.txt", 'w') #open the output file
outfile.write(findAssembly(nodes_list)) #write out the file
outfile.close()
#print(findAssembly(nodes_list))
    
