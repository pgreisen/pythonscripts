from itertools import *
import os
import sys


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
        self.list_of_cst = []
        self.path = './'
        # number of combinations
        self.combi = 4

        self.python_path = "~/pythonscripts/dev/GenerateMultiVariants/"

        self.exe_hc = "python "+self.python_path+"/generate_fasta_tmp.py -f ori_HC.fasta -p "

        self.exe_hc_2 = "python "+self.python_path+"/generate_fasta_tmp.py -f tmp_HC.fasta -p "

        self.modulus = 21000



    def combinations(self,iterable, r):
        # combinations('ABCD', 2) --> AB AC AD BC BD CD
        # combinations(range(4), 3) --> 012 013 023 123
        pool = tuple(iterable)
        n = len(pool)
        if r > n:
            return
        indices = range(r)
        yield tuple(pool[i] for i in indices)
        while True:
            for i in reversed(range(r)):
                if indices[i] != i + n - r:
                    break
            else:
                return
            indices[i] += 1
            for j in range(i+1, r):
                indices[j] = indices[j-1] + 1
            yield tuple(pool[i] for i in indices)


    def setdata(self,datafile):
        with open(datafile) as f:
            for line in f:
                self.list_of_cst.append(line.strip())

    def getUniquePositions(self, list_w_mutations):

        unique = True

        for i in list_w_mutations:
            tmp = i.strip() #split('_')
            # chain = tmp[1]
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
        #
        hc_dummy = 1
        name = ""

        hc_chain = 0

        #NativePositionMutation_chain (S32R_HC )
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


    def main(self):

        inputfile = sys.argv[1]

        combi = sys.argv[2]

        self.combi = int(combi)

        self.setdata(inputfile)

        b = self.combinations(self.list_of_cst  , self.combi)

        cstfiles = []

        
        for j in b:

            # print list(j)
            eval = self.getUniquePositions(list(j))
            print j
            if( eval == True ):
                cstfiles.append(j)

        exe = []
        dummy = 1
        for mutant in cstfiles:
            exe.append( self.generate_fasta(mutant,str(dummy)) )
            dummy += 1

        i = 0
        for line in exe:
            if(i % self.modulus == 0):
                f = open("gen_"+str(i)+".sh",'w')
            f.write(line)

            i += 1
        print "DONE"



if __name__ == "__main__":
   run = Combinations()
   run.main()

