from numpy import mean,sqrt,var
import os,shutil,sys,argparse
import operator
from collections import OrderedDict
from collections import defaultdict
from operator import itemgetter
from numpy import *
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
        self.residues_to_analyse = ["PHE","TYR","TRP","MET","VAL","HIS","ARG","LYS","CYS","PRO","GLN","ASN","GLU","ASP","SER","THR","ILE","LEU"]
        self.rb_cutoff_probability = -0.05


    def set_rb_table(self,file):
        # dummy boolean
        db = False
        flg = 1

        with open( file ,'r') as f:
            for line in f:
                # print flg
                if(len(line) > 19 and line[0:19] == "core.init: command:"):
                    tmp_s = line.split()
                    string_pdb = [s for s in tmp_s if ".pdb" in s]
                    self.pdbfile = string_pdb[0].split('/')[-1]


                if (len(line) > 56 and line[0:56] == self.aa_string):
                    print "Starting analysis..."
                    flg = flg + 1
                    db = True

                if( len(line) > 31 and line[0:31] == self.break_string ):
                    print "Finished with analysis"
                    db = False
                    print "#####################"
                    break


                if( db == True and flg == 2):
                    if(line[0:3] in self.residues_to_analyse ):
                        tmp = line.split()
                        key = tmp[0][0:3]
                        if( key not in self.rbtable.keys() ):

                            self.rbtable[ key ] = []

                        self.rbtable[ key ].append( float ( tmp[2] ) )


    def write_table(self):

        dummy_name = ""
        prbs = []
        tmp_string = ""
        with open("stats.dat",'w') as f:
            f.write("AA,mean,geo-mean,sd,N,\n")
            for key in self.rbtable:
                # print key, self.rbtable[key],len(self.rbtable[key])
                f.write(str(key)+","+str(round( mean(self.rbtable[key]),3))+","+str(round(stats.mstats.gmean( self.rbtable[key]),3))+","+ str(round(sqrt(var(self.rbtable[key])),3))+","+str(len(self.rbtable[key]))+"\n")

            # prbs.append(float(self.rbtable[key]))
            # prbs.append( self.rbtable[key])
            # print self.rbtable[key] , mean( self.rbtable[key] )
            # f.write(key+"\t"+self.rbtable[key]+"\n")
            # if(float(self.rbtable[key]) < self.rb_cutoff_probability ):
            # dummy_name = dummy_name+"FAILED"


        # geom = round(stats.mstats.gmean( asarray( prbs )),3)
        for key in self.rbtable.keys():
            with open("rb_"+key+".dat",'w') as f:
                for line in self.rbtable[key]:
                    f.write(str(line)+'\n')






    def main(self):
        parser = argparse.ArgumentParser(description="Script to analyse output from Rotamer Boltzmann Analysis")
        parser.add_argument('-f',dest='file', help='Data file')
        parser.add_argument('--cutoff',dest='threshold', help='Energy differences below this threshold are considered noise ( default 1 REU )')
        parser.add_argument('--rb_prob',dest='rb_cutoff_probability', help='The rotamer probability cutoff - default is 0.05 ',type=float)

        args_dict = vars( parser.parse_args() )

        for item in args_dict:
            setattr(self, item, args_dict[item])


        # loop over all the files in the directory
        path = './'
        files = os.listdir(path)
        for fl in files:
            if( os.path.isfile(fl) ):
                self.set_rb_table( fl )

        self.write_table()



if __name__ == "__main__":
    run = RBAnalysis()
    run.main()

