import sys

datafile = sys.argv[1]

sequences = {}

with open(datafile,'r') as f:
    for line in f:

        if(line[0] == ">"):
            tmp = line.split()
            ids = tmp[0]
            sequences[ids] = ""
        else:
            tmp = line.strip()
            sequences[ids] = sequences[ids] + tmp
   
for key in sequences:
    with open(key[1:]+'.fasta','w') as f:
        f.write(key+'\n')
        f.write(sequences[key]+'\n')
        
