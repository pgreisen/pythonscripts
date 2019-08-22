import sys

inputfile = sys.argv[1]
muts = []
with open(inputfile, 'r') as f:
    for line in f:
        muts.append(line.strip() )

with open("muts_single.txt",'w') as f:
    for line in set(muts):
        f.write(line+"\n")

