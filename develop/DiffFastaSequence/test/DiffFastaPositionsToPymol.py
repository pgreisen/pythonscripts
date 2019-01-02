#!/usr/bin/env python
import os,shutil, sys, math
from numpy import mean
from pylab import *
import argparse,csv
import operator


class DiffFastaPositionsToPymol:

    def __init__(self):
        self.offset = 0
        self.variantfile = ""
        self.wildtypefile = ""
        self.xtal_offset = 0


    def get_fasta_sequences(self,fastafile):
        seq = {}
        with open(fastafile, 'r') as f:
            for line in f:
                if (line[0] == '>'):
                    key = line.strip()
                    seq[key] = ""
                else:
                    seq[key] = seq[key] + line.strip()

        return seq

    def get_difference_to_wildtype(self,variant,wildtype):
        '''
        :param variant:
        :param wildtype:
        :return:
        '''

        mutations_dict = {}
        pymol_viz_dict = {}
        for i in wildtype.keys():
            seqlength = len(wildtype[i])
            aseq = wildtype[i]
        for var in variant.keys():

            dummy = 1 + self.offset
            mutations = []
            pymol_viz = []
            gap_counter = 0

            bseq = variant[var]

            for i in range(seqlength):
                if (aseq[i] != bseq[i]):
                    mutations.append(aseq[i] + str(dummy) + bseq[i])
                    if (bseq[i] != '-'):
                        pymol_viz.append(str(dummy - gap_counter))
                    else:
                        gap_counter += 1
                dummy += 1
            mutations_dict[var] = mutations
            pymol_viz_dict[var] = pymol_viz

        return mutations_dict,pymol_viz_dict


    def write_pymol_file(self,pymol_dictionary):
        with open("pymol_viz.pml",'w') as f:
            for i in pymol_dictionary.keys():
                tmpname = ""
                tmpresnr = ""
                for j in pymol_dictionary[i]:
                    tmpname += j+"_"
                    tmpresnr += str(j)+"+"

                f.write("create "+i[1:]+"_"+tmpname[0:-1]+", resi "+tmpresnr[0:-1]+'\n')


    def main(self):


        parser = argparse.ArgumentParser(description="Analysis differences between variants and wildtype sequence")

        parser.add_argument('-v',dest='variantfile', help='Fasta file with multiple variant sequences')
        parser.add_argument('-w', dest='wildtypefile', help='Fasta file with reference sequence')
        parser.add_argument('-n', dest='xtal_offset', help='Offset numbering between crystal structure and fasta sequence')

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        variant = self.get_fasta_sequences( self.variantfile )
        wildtype = self.get_fasta_sequences(self.wildtypefile)
        muts,viz = self.get_difference_to_wildtype( variant, wildtype)

        self.write_pymol_file(viz)




if __name__ == "__main__":
   run = DiffFastaPositionsToPymol()
   run.main()




'''
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
'''