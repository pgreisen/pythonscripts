from itertools import *
import os
import sys
import numpy as np
from multiprocessing import Pool

from concurrent.futures import ThreadPoolExecutor
import threading



from GenerateVariantsinFastafile import GenerateVariantsinFastafile
sys.path.append("/z/insulin/users/pjug/pythonscripts/dev/mAbs_GenerateFastaDesign")

'''

python combine_positions.py INPUTFILE


where INPUTFILE contains the aa, positions, and chains to be combined - the syntax in the file is as
following:

S32R_LC
S32H_HC
etc

Next, set the number of combinations e.g. 2 and the modulus can be set such that it is divided into
multiple files.


'''



class Combinations:

    def __init__(self):
        # contain mutations
        self.list_of_mutations = []
        self.path = './'
        self.sequence_variants = {}
        # number of combinations
        self.combi = 4
        self.debug = 0
        self.python_path = "." # "~/pythonscripts/develop/GenerateMultiVariants"

        self.exe_hc = "python "+self.python_path+"/generate_fasta_tmp.py -f ori_HC.fasta -p "

        self.exe_hc_2 = "python "+self.python_path+"/generate_fasta_tmp.py -f tmp_HC.fasta -p "
        self.init_name = ""
        self.modulus = 21000
        self.fastasequence = ""
        self.processors = 8
        self.gvf = GenerateVariantsinFastafile()

    def combinations(self,iterable, r):
        # combinations('ABCD', 2) --> AB AC AD BC BD CD
        # combinations(range(4), 3) --> 012 013 023 123
        pool = tuple(iterable)
        n = len(pool)
        if r > n:
            return
        indices = list(range(r))
        yield tuple(pool[i] for i in indices)
        while True:
            for i in reversed(list(range(r))):
                if indices[i] != i + n - r:
                    break
            else:
                return
            indices[i] += 1
            # python2
            # for j in range(i+1, r):
            # python3
            for j in list(range(i + 1, r)):
                indices[j] = indices[j - 1] + 1
            yield tuple(pool[i] for i in indices)

    def setdata(self,datafile):
        '''
        Make sure that there are no redundant mutations in the input
        '''
        with open(datafile) as f:
            for line in f:
                self.list_of_mutations.append(line.strip() ) 
        self.list_of_mutations = list(set(self.list_of_mutations ))

    def getUniquePositions(self, list_w_mutations):
        unique = True
        for i in list_w_mutations:
            tmp = i.strip()
            position = tmp[1:-1]
            for j in list_w_mutations:

                if(i != j):
                    tmp2 = j.strip()

                    position2 = tmp2[1:-1]

                    if( position == position2 ):
                        # print "same:   ",chain, chain2, position, position2, list_w_mutations
                        return False
                    else:
                        continue

        return unique


    def generate_fasta(self,fastafile,groupnr):
        exe = ""
        hc_dummy = 1
        name = ""
        hc_chain = 0

        for tmpposition in fastafile:

            tmpstr = tmpposition.strip()
            tmp_native = tmpstr[0]
            tmp_pos = tmpstr[1:-1]
            tmp_mutant = tmpstr[-1]
            tmp_chain = "A"

            if( hc_dummy == 1 ):
                tmp_string = self.exe_hc+tmp_pos+" -c "+\
                         tmp_chain+" -n "+tmp_native+" -m "+tmp_mutant+" -g "+groupnr+"\n mv tmp.fasta tmp_HC.fasta\n"

                hc_chain = 1
                exe = exe  + tmp_string

            if( hc_dummy > 1):

                tmp_string = self.exe_hc_2+tmp_pos+" -c "+\
                         tmp_chain+" -n "+tmp_native+" -m "+tmp_mutant+" -g "+groupnr+"\nmv tmp.fasta tmp_HC.fasta\n "
                exe = exe  + tmp_string
                hc_chain = 1
            hc_dummy += 1
            name = name + tmpposition+"_"

        name_hc = name +"HC_"+groupnr+".fasta"
        tmp_move_hc = "mv tmp_HC.fasta "+name_hc+"\n"
        exe = exe +tmp_move_hc+"echo \"group id: \" "+groupnr+"\n"
        return exe

    def generate_fasta_in_parallel(self, mutationalstr, groupnr):
        exe = ""
        hc_dummy = 1
        name = self.init_name
        hc_chain = 0
        newfastasequence = self.fastasequence
        for tmpposition in mutationalstr:
            tmpstr = tmpposition.strip()
            tmp_native = tmpstr[0]
            tmp_pos = int(tmpstr[1:-1])
            tmp_mutant = tmpstr[-1]
            tmp_chain = "A"
            newfastasequence = self.gvf.get_mutated_fasta_string( tmp_pos, tmp_native, tmp_mutant, newfastasequence)
            name += tmpposition+"_"

        name += str(groupnr)

        return name, newfastasequence

    def get_generator_of_combinations(self):
        b = self.combinations(self.list_of_mutations, self.combi)
        return b


    def setup(self, inputfile, combinations,fastafile):
        """
        :param inputfile:
        :param combinations:
        :return:
        """
        self.combi = int(combinations)
        self.setdata(inputfile)
        b = self.combinations(self.list_of_mutations, self.combi)

        mutations = []
        self.fastasequence = self.gvf.get_fastasequence_from_file( fastafile)
        # insert method here to evaluate the WT AA with the given fasta
        self.evalute_fasta_file()

        for j in b:
            eval = self.getUniquePositions(list(j))
            if(self.debug == 1):
                print(j)
            if (eval == True):
                mutations.append(j)

        # rewrite to parallel
        dummy = list(range(len(mutations)) )
        # print("length of designs:: ",dummy)

        with ThreadPoolExecutor(self.processors) as p:
            name_seq = p.map(self.generate_fasta_in_parallel,mutations,dummy)

        with open("designs_"+str(self.combi)+".fasta",'w') as f:
            for i in name_seq:
                f.write(">"+i[0]+"\n")
                f.write(i[1]+"\n")

    def evalute_fasta_file(self):
        for eval_wt in self.list_of_mutations:
            print(eval_wt)
            # count from 0 in python
            pos = int(eval_wt[1:-1] ) -1
            aa = eval_wt[0]
            assert aa == self.fastasequence[pos], "WT =  "+self.fastasequence[pos]+str(pos)+" input AA: "+aa

    def main(self):

        inputfile = sys.argv[1]
        combi = sys.argv[2]
        fastafile = sys.argv[3]
        self.setup(inputfile,combi,fastafile)


if __name__ == "__main__":
   run = Combinations()
   run.main()

