import re, sys
# import fasta file

fastafile = sys.argv[1]

fasta_seqs = {}

with open(fastafile, 'r') as f:
    for line in f:
        if(line[0] == ">"):
            key = line[1:].strip()
            fasta_seqs[key] = ""
        else:
            fasta_seqs[key] += line.strip()

for i in fasta_seqs.keys():
    result = []
    result.append( re.findall(r'N[A-Z][ST]',fasta_seqs[i] ) )
    if(len(result[0]) > 0):
        print("N-Glc present: ",result,i)
