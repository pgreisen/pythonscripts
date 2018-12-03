#!/usr/bin/env python
import sys
import argparse,csv
## Sequence difference between the model

class AnalyseMutations:


    def __init__(self):
        self.native = ""
        self.design = ""
        self.chain = 'A'
    
    # Requires a pdb file
    # Returns sequence 
    def get_seq(self,filename):
        tmp = open(filename,'r')
        resid = ''
        seq = []
        j = 0
        for line in tmp:
            tm = line.split()
            if tm[0] == 'ATOM' and j == 0:
                resid = tm[5]
                seq.append(tm[3]+' '+tm[5])
                j = j+1
            if tm[0] == 'ATOM' and tm[5] != resid:
                seq.append(tm[3]+' '+tm[5])
                resid = tm[5]
        return seq


    def get_seq_pdbfile(self,filename):
        tmp = open(filename,'r')
        resid = '0'
        seq = []
        j = 0
        for line in tmp:
            tm = line.split()
            if line[0:4] == 'ATOM' and resid != line[22:26]:
                resid = line[22:26]
                seq.append(line[17:20]+' '+tm[5])
        return seq

    # Requires two arrays (same length right now)
    # Wild type = a and b = rosetta output
    # Prints the differences between them
    # Prints number of mutations
    # Returns the a list of the mutant's mutations
    def print_def(self,a,b):
        # write pml file
        string = ''
        mt = []
        ln = len(a)
        i = 0
        mut = 0
        tmp_str = ""
        for i in range(ln):
            if a[i] != b[i]:
                mut = mut +1
                # tmp string
                wt_tmp_aa,wt_tmp_pos = a[i].split() 
                wt_tmp_aa = self.get_single_letter_aa( wt_tmp_aa )
                mut_tmp=self.get_single_letter_aa( b[i].split()[0] )
                tmp_str += wt_tmp_aa+wt_tmp_pos+mut_tmp+","
                tmp = a[i].split()
                string = string + '+' + str(tmp[1])


        pymol_file = open("mut.pml",'w')
        pymol_file.write("load "+self.design+", "+self.design[0:-4]+"\n")
        pymol_file.write("show cartoon\n")
        pymol_file.write("hide lines\n")
        pymol_file.write("show sticks, het\n")
        pymol_file.write("color cyan, name C* and het\n")
        pymol_file.write("create "+tmp_str[0:-1].replace(',','_')+", resi "+string+" and chain "+self.chain+" and "+self.design[0:-4]+" \n")
        pymol_file.write("show sticks, "+tmp_str[0:-1].replace(',','_')+"\n")
        pymol_file.write("color copper, "+tmp_str[0:-1].replace(',','_')+" and name C*\n")
        pymol_file.write("hide everything, elem h\n")
        pymol_file.write("center het")
        pymol_file.close()
        return tmp_str


    # get number of mutations
    def get_number_of_mutations(self,a,b):
        ln = len(a)
        assert ln == len(b)
        i = 0
        mut = 0
        for i in range(ln):
            if a[i] != b[i]:
                mut = mut +1
        return mut

    def get_single_letter_aa(self,threeletter):
        res_map = {
            'CYS' : 'C',
            'HIS' : 'H',
            'ASP' : 'D',
            'GLU' : 'E',
            'LYS' : 'K',
            'SER' : 'S',
            'ARG' : 'R',
            'THR' : 'T',
            'ASN' : 'N',
            'GLN' : 'Q',
            'TYR' : 'Y',
            'MET' : 'M',
            'GLY' : 'G',
            'TRP' : 'W',
            'PRO' : 'P',
            'ALA' : 'A',
            'VAL' : 'V',
            'ILE' : 'I',
            'LEU' : 'L',
            'PHE' : 'F'
            }
        return res_map[threeletter]

    def get_mutations(self, filea, fileb):
        seq_a = self.get_seq_pdbfile( filea )
        seq_b = self.get_seq_pdbfile( fileb )
        return self.get_number_of_mutations(seq_a, seq_b)


    def main(self):
        parser = argparse.ArgumentParser(description="Analyse mutations between two pdb files")
        parser.add_argument("-n", "--native", dest="native", help="Mutations relative to this input")
        parser.add_argument("-d", "--design", dest="design", help="Mutations for this input" )
        parser.add_argument("-c", "--chain", dest="chain", default="A",help="chain with the mutations" )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        a = self.get_seq_pdbfile(self.native)
        b = self.get_seq_pdbfile(self.design)

        seq_a = len(a)
        seq_b = len(b)

        #print("Length of "+self.native+" is", seq_a)
        #print("Length of "+self.design+" is", seq_b)
        #print("Seq1 : "+self.design+" Seq2 : "+self.native)
        mutational_string = self.print_def(a,b)
        print(self.design+": "+mutational_string[0:-1])


if __name__ == "__main__":
   analysemutations = AnalyseMutations()
   analysemutations.main()
