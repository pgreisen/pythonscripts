from numpy import mean,sqrt,var
import os,shutil,sys,argparse
import operator
from collections import OrderedDict
from collections import defaultdict

'''

Take an output file from greedy optimization and sorts the data

'''

class AnalyseGreedyTable:

    def __init__(self):
        self.file = ""
        self.threshold = 1
        self.greedytable = defaultdict(list)
        self.average_energy = []


    def get_file(self):
        global lowest_mutations
        lowest_mutations = []
        with open(self.file,'r') as f:
            tmp_native_residue = ""
            for line in f:
                # get residue number
                residue_nr = str(line.split()[1:2]).lstrip('\[\"\[\'\(\'')
                tmp_line = line.split()[3:]
                tmp_energy = []
                native_residue = 0.0
                lowest_energy_residue = 0.0
                tmp_lowest_mutations = ""
                native = False
                native_score = 0.0
                # loop over data file
                for i in tmp_line:
                    tmp_line2 = i.split(':')
                    # identify native residue
                    if(str(tmp_line2[0])[-1] == '*'):
                        native_residue = float(tmp_line2[1])
                        tmp_native_residue = tmp_line2
                        # key for table
                        key = str(tmp_native_residue[0])+str(residue_nr).rstrip('\'\]')
                        # print "key for dictionary: ", key.replace('*',"")
                        native = True
                        native_score = float( tmp_line2[1] )

                        #self.greedytable[key].append((str()))

                    if( float(tmp_line2[1]) < lowest_energy_residue ):
                        lowest_energy_residue = float(tmp_line2[1])
                        tmp_lowest_mutations = tmp_line2
                        self.average_energy.append(float(tmp_line2[1]))
                        print native_residue, lowest_energy_residue
                        diff = native_residue - lowest_energy_residue

                        print diff

                        if(diff > self.threshold and native == False):
                            lowest_mutations.append(str(residue_nr)+' '+str(tmp_lowest_mutations)+' Energy difference: '+str(diff))
                            native = False

        return lowest_mutations

    def main(self):
        parser = argparse.ArgumentParser(description="Script to analyse GreedyOptimization")
        parser.add_argument('-f',dest='file', help='Data file')
        parser.add_argument('--cutoff',dest='threshold', help='Energy differences below this threshold are considered noise ( default 1 REU )')

        args_dict = vars( parser.parse_args() )

        for item in args_dict:
            setattr(self, item, args_dict[item])

        lowest_mutations = self.get_file()



if __name__ == "__main__":
    run = AnalyseGreedyTable()
    run.main()

