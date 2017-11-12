from numpy import mean,sqrt,var
import os,shutil,sys,argparse
import operator
from collections import OrderedDict
from collections import defaultdict
from operator import itemgetter
'''

Take an output file from greedy optimization and sorts the data

'''

class AnalyseGreedyTable:

    def __init__(self):
        self.file = ""
        self.threshold = 1
        self.greedytable = defaultdict(list)
        self.average_energy = []
        self.total_mutations = []


    def get_file(self):
        global lowest_mutations
        lowest_mutations = []
        with open(self.file,'r') as f:
            tmp_native_residue = ""

            for line in f:
                key = ""
                # get residue number
                residue_nr = str(line.split()[1:2]).lstrip('\[\"\[\'\(\'')
                tmp_line = line.split()[3:]
                native_score = 0.0
                substitutions = []
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
                        native_score = float( tmp_line2[1] )

                    else:
                        residue = str(tmp_line2[0])+str(residue_nr).rstrip('\'\]')
                        substitutions.append( [ str(residue) , float(tmp_line2[1])] )
                if( len(line) > 1):
                    for i in substitutions:
                        i[1] = round( native_score - i[1], 3)
                    self.greedytable[ key.replace('*',"") ] =  substitutions
        return 0

    def write_to_file(self,total_mutations):
        with open("mutations.txt",'w') as f:
            f.write("Residue,ddREU")
            for line in total_mutations:
                f.write(str(line[0])+","+str(line[1])+"\n")


    def main(self):
        parser = argparse.ArgumentParser(description="Script to analyse GreedyOptimization")
        parser.add_argument('-f',dest='file', help='Data file')
        parser.add_argument('--cutoff',dest='threshold', help='Energy differences below this threshold are considered noise ( default 1 REU )')

        args_dict = vars( parser.parse_args() )

        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.get_file()

        total_mutations = []
        for key in self.greedytable:
            for value in self.greedytable[key]:
                # contains new amino acid and it rosetta score
                print key,value

                new_value = key+value[0][0:3]

                total_mutations.append( (new_value, value[1]) )

        #print self.greedytable

        #
        total_mutations.sort(key=lambda x: x[1], reverse=True)
        # sorted(total_mutations, key=itemgetter(0) )
        #for i in total_mutations:
        #    print i
        # print total_mutations
        self.write_to_file( total_mutations)


if __name__ == "__main__":
    run = AnalyseGreedyTable()
    run.main()

