from numpy import mean,sqrt,var
import sys
from collections import defaultdict
from collections import OrderedDict

class AnalyseGreedyTable:

    def __init__(self):
        self.aa = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        self.aa_values = []
        # times above mean
        self.factor = 1


    def get_sorted_hashtable(self, hashtable):
        return OrderedDict(sorted(hashtable.items(), key=lambda x: x[1],reverse=True)) #[0:maxvalue])


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

    def get_mean_and_sd(self):
        return mean( self.aa_values ), sqrt(var( self.aa_values) )


    def write_aa_above_mean(self,dictionary_w_values):
        d = defaultdict(list)
        mn, sd = self.get_mean_and_sd()
        for k in dictionary_w_values:
            # the value is a new dictionary
            for kk in  dictionary_w_values[k]:
                if( float( dictionary_w_values[k][kk] ) < mn - self.factor*sd ):
                    residue = kk+" "+k
                    d[residue] = str(round( float(dictionary_w_values[k][kk]), 3))

        d = self.get_sorted_hashtable( d )
        with open("Substitutions_below_mean", 'w') as f:
            for key in d:
                f.write(key+": "+d[key]+"\n")

        with open("generate_substitutions.sh", 'w') as f:
            for key in d:
                resname, resnr = key.split()
                tmp_string = "sh ~greisen/files/MutateResidue/run_residue_number.sh $1 $2 "+str(resnr)+" "+str(resname)
                f.write(tmp_string+"\n")
                f.write("mv *_0001.pdb "+str(resname)+"_"+str(resnr)+".pdb\n")





    def main(self):
        table = sys.argv[1]    
        tmp = open(table,'r')
        whole_dictionary = {}
        for line in tmp:
            if ( len(line)  < 6 ):
                continue
            dictionary = {}
            residue_nr = str(line.split()[1:2]).lstrip('\[\"\[\'\(\'')
            residue_nr = residue_nr.rstrip('\]\"\]\'\]\'')
            tmp_line = line.split()[3:]
            for i in tmp_line:

                if( len(i) < 2):
                    continue

                tmp_line2 = i.split(':')
                amino_acid = str(tmp_line2[0])[0:3]
                value = tmp_line2[1]
                dictionary[amino_acid] =  value
                self.aa_values.append( float( value ) )
            whole_dictionary[residue_nr] = dictionary

        self.write_matrix_to_file( whole_dictionary )
        self.write_aa_above_mean(whole_dictionary)


if __name__ == "__main__":
   run = AnalyseGreedyTable()
   run.main()