import sys, subprocess, os

pdbfile = sys.argv[1]
mutationsfile = sys.argv[2]

# get fastafile
pdb_to_fasta = "perl /Users/pgreisen/bin/pdb_fasta.pl  "+pdbfile+ " > pdb.fasta"
os.system(pdb_to_fasta)
fastaseq = {}

with open("pdb.fasta", 'r') as f:
    for i in f:
        if(i[0] == '>'):
            key = i.strip()
            fastaseq[key] = ""
        else:
            fastaseq[key] += i.strip().replace("-","")
mutations = []
with open(mutationsfile,'r') as f:
    for line in f:
        mutations.append(line.strip())
mutations = list(set(mutations))
fix = []
for i in mutations:
    mutstring = "Mut: "+i[0]+" pos: "+str(i[1:-1])+" WT: "+fastaseq[key][int(i[1:-1])]

    if( i[0] != fastaseq[key][int(i[1:-1])-1 ]):
        fix.append(mutstring)
print(fix)
