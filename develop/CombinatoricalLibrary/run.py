#!/usr/bin/env python
import sys
from numpy import *
import argparse
from random import randint
import pylab as plt
from pylab import *
import pandas as pd
import argparse



class CombinatoricsLibrarySize:

    def __init__(self):
        self.csvfile = ""
        self.pos = 0
        self.AA = 20
        self.samples = 1000000
        self.dummy_sequences = []
        self.native_prb = 60
        self.charged_residues_prb = 7.5
        self.aa_size = 9
        self.nsamples = 1
        self.number_of_native_positions = 2
        self.number_of_charged_residues = 3
        self.charged_bias = self.native_prb + self.charged_residues_prb * self.number_of_charged_residues


    def get_mutation_uniform(self, randomnumber):
        if (randomnumber <= self.native_prb):
            return "N"
        else: # (randomnumber > self.native_prb and randomnumber <= 100):
            return "X"

    def get_mutation(self, randomnumber):

        if (randomnumber <= self.native_prb + self.charged_residues_prb):
            return "N"

        else:  # elif (randomnumber > self.native_prb and randomnumber <= 100):
            return "X"

    def get_mutation_charged_residue_wt(self, randomnumber):
        # position 86,89 prb = 0.9 for WT and charged residue

        if (randomnumber <= self.native_prb + self.charged_residues_prb):
            return "N"

        elif(randomnumber > self.native_prb and randomnumber <= self.native_prb + 3*self.charged_residues_prb):
            return "C"

        else:  # elif (randomnumber > self.native_prb and randomnumber <= 100):
            return "X"


    def get_mutation_charged_residue(self, randomnumber):
        # position 86,89 prb = 0.9 for WT and charged residue

        if (randomnumber <= self.native_prb):
            return "N"

        elif(randomnumber > self.native_prb and randomnumber <= self.charged_bias):
            return "C"

        else:  # elif (randomnumber > self.native_prb and randomnumber <= 100):
            return "X"


    def mutation_analysis(self):

        mutations_analysis = []
        for i in range(self.samples):
            tmp_seq = []
            # make sequences
            for i in range(self.aa_size):

                # generate random number
                rdv = randint(0, 100)

                # positions with charged aa
                if (i <= self.number_of_native_positions):

                    tmp_seq.append(self.get_mutation(rdv))
                else:
                    tmp_seq.append(self.get_mutation_uniform(rdv))

                    # [1, 2, 3, 4, 1, 4, 1].count(1)
            mutations_analysis.append(tmp_seq.count("X"))

        return mutations_analysis

    def mutational_analysis_w_charges(self):

        mutations_w_charge = []
        for i in range(self.samples):
            tmp_seq = []
            # make sequences
            for i in range(self.aa_size):
                # generate random number
                rdv = randint(0, 100)
                # positions with charged aa
                if (i <= self.number_of_native_positions):
                    tmp_seq.append(self.get_mutation_charged_residue_wt(rdv))
                else:
                    tmp_seq.append(self.get_mutation_charged_residue(rdv))
                    # [1, 2, 3, 4, 1, 4, 1].count(1)
            mutations_w_charge.append(tmp_seq.count("C"))

        return mutations_w_charge


    def get_dictionary_w_probabilities(self,df):
        dict_ = {}
        df.fillna(value=0,inplace=True)
        for i,j,k in zip(df.pos,df.wt,df.group):
            dict_[i] = []
            dict_[i].append(float(j))
            dict_[i].append(float(j)+float(k))
            self.pos += 1

        return dict_


    def compute_probabilities(self,dict):
        import random
        # random.uniform(0, 1)
        assert len(dict.keys()) == self.pos
        # for n in range(self.nsamples):
        mutational_freq = []
        for i in range(self.samples):
            muts = ""
            for j in dict.keys() :
                random_number = random.uniform(0, 1)
                if( random_number < dict[j][0] ):
                    muts += "N"
                elif( random_number < (dict[j][1]) ):
                    muts += "G"
                else:
                    muts += "X"
            g_ = muts.count("G")
            x_ = muts.count("X")
            tot_mut = g_ + x_
            mutational_freq.append( tot_mut)
        return mutational_freq



    def nCr(self,n, r):
        import math
        f = math.factorial
        return f(n) / f(r) / f(n - r)


    def write_to_file(self,mutational_analysis):
        with open("library.csv",'w') as f:
            f.write(",Counts,Frequency,Combinations,Combi_frac\n")
            for i in range(self.pos + 1):

                combi = self.nCr(self.pos,i)*20**i
                combi_frac = str(round(  mutational_analysis.count(i) / combi, 4 ) )
                f.write(str(i)+" mutations:,"+str(mutational_analysis.count(i))+","+str(round(float(mutational_analysis.count(i))/self.samples,2))+","+str(combi)+","+combi_frac+"\n" )
            f.write("Mean: " + str( mean(mutational_analysis) )+ "\n")
            f.write("Mean: " + str( median(mutational_analysis) )+ "\n")


    def main(self):
        mutational_analysis = self.mutation_analysis()
        parser = argparse.ArgumentParser(description=" ")
        parser.add_argument("-f", "--file", dest="csvfile", help="CSV file with probabilities", default="")
        parser.add_argument("-a", "--aminoacids", dest="AA", help="Number of amino acids available in library", default=20,type=int)
        parser.add_argument("-s", "--samplesize", dest="samples", help="How big a population do we sample", default=1000000000,type=int)
        parser.add_argument("-n", "--numberofsamples", dest="nsamples", help="How many samples should we simulate", default=1,type=int)
        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        df = pd.read_csv(self.csvfile)
        dict = self.get_dictionary_w_probabilities(df)

        mutational_analysis = self.compute_probabilities(dict)
        print mean(mutational_analysis)
        print median(mutational_analysis)
        hist(mutational_analysis)
        savefig("mutations.png")
        self.write_to_file(mutational_analysis)

        #mutations_w_charge = self.mutational_analysis_w_charges()
        #print mean(mutations_w_charge)
        ##print median(mutations_w_charge)
        #clf()
        ##hist(mutations_w_charge)
        #savefig("mutations_w_charge.png")
        #
        #with open("charge_library.csv", 'w') as f:
        #    for i in range(self.aa_size):
        #        f.write(str(i) + " mutations:," + str(mutations_w_charge.count(i)) + "," + str(round(float(mutations_w_charge.count(i)) / self.samples, 2)) + "\n")

if __name__ == "__main__":
    run = CombinatoricsLibrarySize()
    run.main()
