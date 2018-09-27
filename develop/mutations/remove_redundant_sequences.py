#!/usr/bin/env python
import sys,os,shutil
## Sequence difference between the model


# Requires a pdb file
# Returns sequence 
def get_seq(filename):
    tmp = open(filename,'r')
    resid = ''
    seq = []
    j = 0
    for line in tmp:
        tm = line.split()
        if tm[0] == 'ATOM' and j == 0:
            resid = tm[5]
            seq.append(tm[3]+' '+tm[5])
        j = j+1
        if tm[0] == 'ATOM' and tm[5] != resid:
            seq.append(tm[3]+' '+tm[5])
            resid = tm[5]
    return seq


def get_seq_pdbfile(filename):
    tmp = open(filename,'r')
    resid = '0'
    seq = []
    j = 0
    for line in tmp:
        tm = line.split()
        if line[0:4] == 'ATOM' and resid != line[22:26]:

            resid = line[22:26]
            seq.append(line[17:20]+' '+tm[5])

    return seq



# Requires two arrays (same length right now)
# Wild type = a and b = rosetta output
# Prints the differences between them
# Prints number of mutations
# Returns the a list of the mutant's mutations
def print_def(a,b):
    string = ''
    mt = []
    ln = len(a)
    i = 0
    mut = 0

    if( len(a) != len(b) ):
        mut = -1

    else:

        for i in range(ln):

            if a[i] != b[i]:

                mut = mut +1

                mt.append(b[i])
                
                tmp = a[i].split()

                string = string + '+' + str(tmp[1])

    return mut


def main():
    path = './'

    files = os.listdir(path)

    # Make directory for redundant structure
    dst = "Redundant_sequences"

    try:
        os.mkdir(dst)

    except:
        print("Directory already exists")
    
    for fl1 in files:

        if( os.path.isfile(fl1) and fl1[-3:] == 'pdb'):

            a = get_seq_pdbfile(fl1)

            for fl2 in files:

                if(fl1 == fl2):
                    continue

                elif( os.path.isdir(fl2) ):
                    continue

                elif( os.path.isfile(fl2) and fl2[-3:] == 'pdb' ):

                    b = get_seq_pdbfile(fl2)
                    
                    # Contains number of mutations

                    mut = print_def(a,b)

                    if(mut == 0):
                        print("The files have same sequence ", fl1, fl2)
                        shutil.move(fl2,dst)
                else:
                    continue


if __name__ == "__main__":
    main() 
