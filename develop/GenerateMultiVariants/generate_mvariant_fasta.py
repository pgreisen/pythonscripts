#!/usr/bin/env python
import os,shutil, commands, sys, math
from numpy import mean
from pylab import *
import argparse,csv
import operator

'''

Assume a single data file first column with names and second column values:

python generate_fasta.py -f wt.fasta -n NATIVE -p POSITION -m MUTATION -c A -g GROUP_ID


'''

class PlotScorefileRosetta:


    def __init__(self):
        self.data = {}
        self.fastaseq = ""
        self.newfasta = ""
        self.chain = "A"
        self.group = ""
        self.name = ""
        self.filename_id = ""
        self.fasta_header = ""


    def getdata(self, inputfile):
        fasta_seq = []
        with open(inputfile,'r') as f:
            for line in f:
                if( line[0] == '>' ):
                    continue
                else:
                    self.fastaseq = self.fastaseq + line.strip()

    def insert_mutations(self, position, native, newmutation, newfasta ):
        for aa in range (len(self.fastaseq)):
            if aa == position-1:
                assert self.fastaseq[aa] == native

                newfasta = newfasta[0:position-1]+newmutation+newfasta[position:]

        return newfasta


    def write2file(self,filename):
        filename = filename.replace("__", "_")
        self.fasta_header = self.fasta_header.replace("__", "_")
        # print filename
        print self.name
        if(self.name != "-1"):
            filename = self.fasta_header+".fasta"

        with open(filename[0:]+self.filename_id,'w') as f:
            f.write(">"+self.fasta_header+'\n')
            f.write(self.newfasta)


    def main(self):

        '''


        python generate_fasta.py -f wt.fasta -n NATIVE -p POSITION -m MUTATION -c LC -g GROUP_ID


        '''



        parser = argparse.ArgumentParser(description="Insert mutations into fasta file ")

        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="scorefile", help="Initial Rosetta score file")

        parser.add_argument("-m", "--mutations", dest="mutation", help="New amino acid entity")

        parser.add_argument("-n", "--native", dest="native", help="Old amino acid entity")
      
        parser.add_argument("-p", "--position", dest="position", help="Position to mutate" )

        parser.add_argument("-c", "--chain", dest="chain", help="Light or heavy chain",default="" )

        parser.add_argument("-g", "--group", dest="group", help="Group to determine heavy and light chain pairing", default="")

        parser.add_argument("-x", "--name", dest="name", help="Name for fasta file and header", default="-1")

        parser.add_argument("-y", "--id", dest="id", help="Name for fasta file and header", default="")


        input_variables = parser.parse_args()

        self.filename_id = input_variables.id

        self.group = input_variables.group


        self.chain = input_variables.chain

        self.getdata( input_variables.scorefile )

        self.newfasta = self.fastaseq

        self.name = input_variables.name

        tmp_filename = ""

        mutations = input_variables.mutation.split(",")
        positions = input_variables.position.split(",")
        natives = input_variables.native.split(",")

        for i,j,k in zip(positions,natives, mutations):
            self.newfasta = self.insert_mutations( int(i), j, k, self.newfasta )
            tmp_filename = tmp_filename+"_"+j+str(i)+k+"_"
        filename = input_variables.scorefile.split('.')[0].split('_')[0]+tmp_filename+input_variables.chain+"_"+self.group+".fasta"

        self.fasta_header = tmp_filename[1:-1]


        self.write2file( filename )


if __name__ == "__main__":
   run = PlotScorefileRosetta()
   run.main()
