import sys


def getSeq(datafile, seq):
    with open(datafile,'r') as f:
        for line in f:
            if(line[0] == '>'):
                key = line.split()[0]
                seq[key] = ""
            else:
                seq[key] = seq[key] + line.strip()
    return seq, key


seq = {}

seq1 = sys.argv[1]
seq2 = sys.argv[2]


seq,mutkey = getSeq(seq1,seq)
seq,nativekey = getSeq(seq2,seq)


tmp_seq = []
for key in seq:
    tmp_seq.append(seq[key])

aseq = tmp_seq[0]
bseq = tmp_seq[1]

dummy = 1

seqlength = len(aseq)
mutations = []

for i in range( seqlength ):

    if(aseq[i] != bseq[i]):
        mutations.append(aseq[i]+str(dummy)+bseq[i])
    dummy += 1

mut = ""
for j in mutations:
    mut = mut +","+ j
print(mutkey,mut)
