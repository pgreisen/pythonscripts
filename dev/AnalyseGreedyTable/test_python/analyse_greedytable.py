from numpy import mean,sqrt,var
import os,shutil,sys,argparse
import operator
from collections import OrderedDict

'''

Take an output file from greedy optimization and sorts the data

'''

class AnalyseGreedyTable:

    def __init__(self):
        self.file = ""


    def get_file(self):
        
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
                
                for i in tmp_line:
                    tmp_line2 = i.split(':')

                    # if native residue
                    if (str(tmp_line2[0])[-1] == '*'):
                        native_residue = float(tmp_line2[1])
                        tmp_native_residue = tmp_line2

                    if float(tmp_line2[1]) < lowest_energy_residue:
                        lowest_energy_residue = float(tmp_line2[1])
                        tmp_lowest_mutations = tmp_line2
                        average_energy.append(float(tmp_line2[1]))
                        diff = native_residue - lowest_energy_residue


                        if(diff > 1):
                            # print diff, native_residue, lowest_energy_residue
                            lowest_mutations.append(str(residue_nr)+' '+str(tmp_lowest_mutations)+' Energy difference: '+str(diff))
                            lw_mutations.write(str(tmp_native_residue[0])+' '+str(residue_nr).rstrip('\'\]')+' '+str(tmp_lowest_mutations[0])+' Energy difference: '+str(diff)+'\n')


    def main(self):
        parser = argparse.ArgumentParser(description="Script to analyse GreedyOptimization")
        parser.add_argument('-f',dest='file', help='Data file')

        args_dict = vars( parser.parse_args() )

        for item in args_dict:
            setattr(self, item, args_dict[item])



        lw_mutations = open('lowest_mutations','w')
        

        
        average_energy = []

        lowest_mutations = []


print 'Mean of sample: ',mean(average_energy)
print 'SD of sample: ',sqrt(var(average_energy))
