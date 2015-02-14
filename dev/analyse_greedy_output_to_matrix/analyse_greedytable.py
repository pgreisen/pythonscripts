from numpy import mean,sqrt,var
import sys

class AnalyseGreedyTable:

    def __init__(self):
        self.aa = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

    def write_matrix_to_file(self,datamatrix):
        matrix = open("data.csv",'w')
        matrix.write("Position,")
        for i in self.aa:
            # to add a new line
            if( i != "VAL" ):
                matrix.write(i+',')
            else:
                matrix.write(i+"\n")
        for key, value in datamatrix.iteritems():
            if(key == ''):
                continue
            matrix.write(key+',')
            for aa in self.aa:
                if( value.has_key( aa ) ):
                    matrix.write(value[aa]+',')
                else:
                    matrix.write('NaN,')


    def main(self):
        lw_mutations = open('lowest_mutations','w')    
        table = sys.argv[1]    
        tmp = open(table,'r')
        whole_dictionary = {}
        for line in tmp:
            dictionary = {}
            residue_nr = str(line.split()[1:2]).lstrip('\[\"\[\'\(\'')
            residue_nr = residue_nr.rstrip('\]\"\]\'\]\'')
            tmp_line = line.split()[3:]
            for i in tmp_line:
                tmp_line2 = i.split(':')
                amino_acid = str(tmp_line2[0])[0:3]
                value = tmp_line2[1]
                dictionary[amino_acid] =  value 
            whole_dictionary[residue_nr] = dictionary
        self.write_matrix_to_file( whole_dictionary )


if __name__ == "__main__":
   run = AnalyseGreedyTable()
   run.main()
