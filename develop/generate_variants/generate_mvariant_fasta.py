#!/usr/bin/env python
import os,shutil,  sys, math #, commands
from numpy import mean
from pylab import *
import argparse,csv
import operator
from multiprocessing import Pool

'''

Assume a single data file first column with names and second column values:


python plot_data_w_mean_sd.py -f FILE --histogram True -n Histogram

'''

class GenerateVariantsinFastafile:


    def __init__(self):
        self.data = {}
        self.fastaseq = ""
        self.newfasta = ""
        self.chain = "LC"
        self.group = ""
        self.name = "-1"
        self.filename_id = ""

    def getdata(self, inputfile):
        fasta_seq = []
        with open(inputfile, 'r') as f:
            for line in f:
                if (line[0] == '>'):
                    continue
                else:
                    self.fastaseq = self.fastaseq + line.strip()

    def get_fastasequence_from_file(self, inputfile):
        fasta_seq = ""
        with open(inputfile, 'r') as f:
            for line in f:
                if (line[0] == '>'):
                    continue
                else:
                    fasta_seq += line.strip()
        return fasta_seq

    def insert_mutations(self, position, native, newmutation, newfasta ):
        '''
        :param position:
        :param native:
        :param newmutation:
        :param newfasta:
        :return: fastafile with mutation inserted
        '''
        # print "Input to method: ", position, native, newmutation, newfasta
        for aa in range (len(self.fastaseq)):
            # print aa, position-1
            if aa == position-1:
                # print self.fastaseq[0:position-1]
                # print newmutation
                # print self.fastaseq[position:]
                assert self.fastaseq[aa] == native
                newfasta = newfasta[0:position-1]+newmutation+newfasta[position:]

        return newfasta

    def get_mutated_fasta_string(self, position, native, newmutation, fastasequence):

        '''
        :param position:
        :param native:
        :param newmutation:
        :param newfasta:
        :return: fastafile with mutation inserted
        '''
        newfasta = fastasequence
        for aa in range(len(fastasequence)):
            if aa == position-1:
                # print("Native,Pos,Design",native,position,fastasequence[aa])
                assert fastasequence[aa] == native
                newfasta = newfasta[0:position-1]+newmutation+newfasta[position:]

        return newfasta



    def write2file(self,filename):
        filename = filename.replace("__","_")
        # print filename
        # print self.name
        if(self.name != "-1"):
            filename = self.name+self.group+".fasta"
            # print filename

        with open(filename+self.filename_id,'w') as f:
            f.write(">"+filename.split('.')[0]+'\n')
            f.write(self.newfasta)


    def run_analysis(self,fastafile,mutation, native, position,chain):
        '''
        :param fastafile: native fastafile - i think need to be a sequnce
        :param mutations: mutations string separated with comma
        :param native: native amino acids
        :param position: positions to mutate
        :param chain: heavy or light chain
        :return: fasta sequence with new mutations
        '''
        # self.fastaseq
        self.getdata( fastafile )


        mutations = mutation.split(",")
        positions = position.split(",")
        natives = native.split(",")

        for i, j, k in zip(positions, natives, mutations):
            # print "Debug from method 2: ", i,j,k
            self.fastaseq = self.insert_mutations( int(i), j, k, self.fastaseq)
        return self.fastaseq



    def main(self):

        '''


        python generate_fasta.py -f 1886_lc.fasta -n NATIVE -p POSITION -m MUTATION -c LC -g GROUP_ID


        '''



        parser = argparse.ArgumentParser(description="Insert mutations into fasta file ")

        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="scorefile", help="Initial Rosetta score file")

        parser.add_argument("-m", "--mutations", dest="mutation", help="New amino acid entity")

        parser.add_argument("-n", "--native", dest="native", help="Old amino acid entity")

        parser.add_argument("-p", "--position", dest="position", help="Position to mutate" )

        parser.add_argument("-c", "--chain", dest="chain", help="Light or heavy chain" )

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
        self.write2file( filename )


if __name__ == "__main__":
   run = GenerateVariantsinFastafile()
   run.main()
