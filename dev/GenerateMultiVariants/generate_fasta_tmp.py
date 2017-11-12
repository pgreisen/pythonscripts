#!/usr/bin/env python
import os,shutil, commands, sys, math
from numpy import mean
# from pylab import *
import argparse,csv
import operator

'''

Assume a single data file first column with names and second column values:


python plot_data_w_mean_sd.py -f FILE --histogram True -n Histogram

'''

class PlotScorefileRosetta:


    def __init__(self):
        self.data = {}
        self.fastaseq = ""
        self.newfasta = ""
        self.chain = "LC"
        self.group = ""

    def getdata(self, inputfile):
        fasta_seq = []
        with open(inputfile,'r') as f:
            for line in f:
                if( line[0] == '>' ):
                    continue
                else:
                    self.fastaseq = self.fastaseq + line.strip()
        # print self.fastaseq



    def insert_mutations(self, position, native, newmutation ):
        for aa in range (len(self.fastaseq)):
            if aa == position-1:

                # print self.fastaseq[aa] , native , aa,position
                assert_output = "this fails "+self.fastaseq[aa]+native+str(aa)+"_"+str(position)
                assert self.fastaseq[aa] == native, assert_output

                # print self.fastaseq[0:position-1]
                # print newmutation
                # print self.fastaseq[position:]

                self.newfasta = self.fastaseq[0:position-1]+newmutation+self.fastaseq[position:]


    def write2file(self,filename):
        with open("tmp.fasta",'w') as f:
            f.write(">"+filename.split('.')[0]+'\n')
            f.write(self.newfasta)


    def main(self):

        parser = argparse.ArgumentParser(description="Insert mutations into fasta file ")

        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="scorefile", help="Initial Rosetta score file")

        parser.add_argument("-m", "--mutations", dest="mutation", help="New amino acid entity")

        parser.add_argument("-n", "--native", dest="native", help="Old amino acid entity")
      
        parser.add_argument("-p", "--position", dest="position", help="Position to mutate", type=int )

        parser.add_argument("-c", "--chain", dest="chain", help="Light or heavy chain" )

        parser.add_argument("-g", "--group", dest="group", help="Group to determine heavy and light chain pairing", default="")

        input_variables = parser.parse_args()

        self.group = input_variables.group

        self.chain = input_variables.chain

        self.getdata( input_variables.scorefile )

        if(input_variables.position > 0 ):

            self.insert_mutations( input_variables.position, input_variables.native, input_variables.mutation[0]   )

            filename = input_variables.scorefile.split('.')[0].split('_')[0]+'_'+input_variables.native+str(input_variables.position)+input_variables.mutation[0]+'_'+input_variables.chain+"_"+self.group+".fasta"

        else:

            self.newfasta = self.fastaseq 

            filename = input_variables.scorefile.split('.')[0].split('_')[0]+'_'+input_variables.chain+"_"+self.group+".fasta"


        self.write2file( filename )

        # print self.newfasta



if __name__ == "__main__":
   run = PlotScorefileRosetta()
   run.main()
