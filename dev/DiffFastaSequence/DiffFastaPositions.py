import sys

seq = {}

offset = 0


inputfile = sys.argv[1]

with open(inputfile,'r') as f:
    for line in f:
        if(line[0] == '>'):
            key = line.split()[0]
            seq[key] = ""
        else:
            seq[key] = seq[key] + line.split()[0]


aseq = []
bseq = []
i = 1
for key in seq:
    if(i == 1):
        aseq = seq[key]
    else:
        bseq = seq[key]
    i += 1

seqlength = len(aseq)
dummy = 1+offset
mutations = []
pymol_viz = []

gap_counter = 0

for i in range( seqlength ):
    if(aseq[i] != bseq[i]):
        mutations.append(aseq[i]+str(dummy)+bseq[i])
        if( bseq[i] != '-' ):
            pymol_viz.append(str(dummy - gap_counter))
        else:
            gap_counter += 1

    dummy += 1

with open("viz.pymol",'w') as f:
    f.write("create diff, resi ")
    for i in pymol_viz:
        f.write(i+'+')
    f.write('\n')
    f.write("show sticks, diff")

with open("diff.txt",'w') as f:
    for i in mutations:
        f.write(i+"\n")
