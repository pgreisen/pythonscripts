'''
Uses header to insert codons into DNA string - remember to rename accordingly.


'''

import sys
sys.path.append("/Users/pgreisen/pythonscripts/develop/insert_mutation_into_DNA_string")

import InsertSubstitutions
ris = InsertSubstitutions.InsertSubstitutions()

dnafile = sys.argv[1]
fastafile = sys.argv[2]
codonFILE = sys.argv[3]


sequences = {}

debug_sequence = []

with open(fastafile, 'r') as f:
    for line in f:
        if (line[0] == ">"):
            tmp = line.split()
            ids = tmp[0]
            sequences[ids] = ""
            debug_sequence.append(ids)
        else:
            tmp = line.strip()
            sequences[ids] = sequences[ids] + tmp


print(len(debug_sequence), len(set(debug_sequence) ) )


import collections
print( [item for item, count in collections.Counter(debug_sequence).items() if count > 1])

dna_seqs = {}

for key in sequences.keys():
    stop_of_list = -1
    mut_header = ""
    tmpheader = key[1:].split('_')
    integer_ = isinstance(tmpheader[-1], int)
    if(integer_ == False):
        stop_of_list = len(tmpheader)
    for i in tmpheader[0:stop_of_list]:
        mut_header += i+"_"
    fastaheader, sequence = ris.run( dnafile, codonFILE, mut_header[0:-1])
    dna_seqs[fastaheader] = sequence

from datetime import date
today = date.today()
# write fasta sequence to file
with open(str(today)+"_aa_sequences.fasta",'w') as f:
    for key,value in dna_seqs.items():
        f.write(">"+key+"\n")
        f.write(value+"\n")
