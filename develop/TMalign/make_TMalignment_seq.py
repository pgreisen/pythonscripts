import sys, shutil, os, subprocess, operator
from collections import defaultdict


'''

for w in sorted(d, key=d.get, reverse=True):
  print w, d[w]


'''

structure_sequences = {}

path = './'
score = defaultdict(list)
files = os.listdir(path)
threshold = 0.5

threshold_query_seq = 0.75
threshold_target_seq = 0.75

for fl in files:
    if fl.endswith('logfile'):
        pdbid = fl.split('_')[0]
        ids = pdbid[0:4]
        # read file and get score
        with open(fl,'r') as f:
            seq_access = False
            dummy = 0
            for line in f:

                if("Name of Chain_1:" in line):
                    key_native = line.replace(".//", "").split()[-1]
                elif("Name of Chain_2" in line):
                    key_target = line.replace(".//", "").split()[-1]
                elif("TM-score=" in line and "normalized by length of Chain_1" in line):
                    if(float( line.split()[1]) < threshold_query_seq):
                        break
                elif ("TM-score=" in line and "normalized by length of Chain_2" in line):
                    if (float(line.split()[1]) < threshold_target_seq):
                        break
                elif("(\":\" denotes aligned residue pairs of d < 5.0 A, \".\" denotes other aligned residues)" in line):
                    seq_access = True

                elif(seq_access == True and line[0] !=':' and line[0] != " " and line[0] !='.'):
                    print line
                    if(dummy == 0):
                        if(key_native not in structure_sequences.keys()):
                            structure_sequences[key_native] = line.strip()
                        dummy += 1
                    elif(dummy == 1):
                         structure_sequences[key_target] = line.strip()
                         dummy += 1
                    else:
                        continue
                else:
                    continue

with open("TMalign_seq.fasta",'w') as f:
    for key in structure_sequences.keys():
        f.write(">"+key+"\n")
        f.write(structure_sequences[key]+"\n")
