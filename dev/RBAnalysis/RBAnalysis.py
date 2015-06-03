from numpy import mean,sqrt,var
import os,shutil,sys,argparse
import operator
from collections import OrderedDict
from collections import defaultdict
from operator import itemgetter

from scipy import stats

'''



'''

class RBAnalysis:

    def __init__(self):
        self.file = ""
        self.threshold = 1
        self.rbtable = defaultdict(list)
        self.average_energy = []
        self.total_mutations = []
        self.pdbfile = ""
        self.aa_string = "protocols.simple_filters.RotamerBoltzmannWeight: Residue"
        self.break_string ="ddG before, after modification:"
        self.residues_to_analyse = ["PHE","TYR","TRP"]
        self.rb_cutoff_probability = 0.05


    def set_rb_table(self):
        # dummy boolean
        db = False
        with open(self.file,'r') as f:
            for line in f:

                if(len(line) > 19 and line[0:19] == "core.init: command:"):
                    tmp_s = line.split()
                    string_pdb = [s for s in tmp_s if ".pdb" in s]
                    self.pdbfile = string_pdb[0].split('/')[-1]


                if( len(line) > 31 and line[0:31] == self.break_string ):
                    print "Finished with analysis"
                    db = False
                    break

                if( len(line) > 56 and line[0:56] == self.aa_string ):
                    print "Starting analysis..."
                    db = True

                if(db == True):
                    if(line[0:3] in self.residues_to_analyse ):
                        tmp = line.split()
                        self.rbtable[tmp[0]] = tmp[2]


    def write_table(self):

        dummy_name = ""
        prbs = []
        tmp_string = ""
        for key in self.rbtable:
            tmp_string = tmp_string+key+"\t"+self.rbtable[key]+"\n"
            prbs.append(float(self.rbtable[key]))
            # f.write(key+"\t"+self.rbtable[key]+"\n")
            if(float(self.rbtable[key]) < self.rb_cutoff_probability ):

                dummy_name = dummy_name+"FAILED"
        geom = round(stats.mstats.gmean(prbs),3)

        with open("rb_"+self.pdbfile+dummy_name,'w') as f:
            for line in tmp_string:
                f.write(line)
            f.write("Geometric Mean is: "+str(geom))





    def main(self):
        parser = argparse.ArgumentParser(description="Script to analyse output from Rotamer Boltzmann Analysis")
        parser.add_argument('-f',dest='file', help='Data file')
        parser.add_argument('--cutoff',dest='threshold', help='Energy differences below this threshold are considered noise ( default 1 REU )')
        parser.add_argument('--rb_prob',dest='', help='The rotamer probability cutrb_cutoff_probabilityoff - default is 0.05 ',type=float)


        args_dict = vars( parser.parse_args() )

        for item in args_dict:
            setattr(self, item, args_dict[item])
        self.set_rb_table()
        self.write_table()



if __name__ == "__main__":
    run = RBAnalysis()
    run.main()

