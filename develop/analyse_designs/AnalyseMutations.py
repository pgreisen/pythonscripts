#!/usr/bin/env python
import sys
## Sequence difference between the model

class AnalyseMutations:
    

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
        pymol_file = open("mut.pml",'w')
        pymol_file.write("hide sticks, not het\n")
        string = ''
        mt = []
        ln = len(a)
        i = 0
        mut = 0

        for i in range(ln):
            if a[i] != b[i]:
                print a[i]+'                  '+b[i]
                mut = mut +1
                mt.append(b[i])

                tmp = a[i].split()
                string = string + '+' + str(tmp[1])

        print 'Number of mutations',mut
        print string
        pymol_file.write("show sticks, resi "+string+"\n")
        pymol_file.write("hide everything, elem h")
        pymol_file.close()
        return mt


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



    def get_mutations(self, filea, fileb):
        seq_a = self.get_seq_pdbfile( filea )
        seq_b = self.get_seq_pdbfile( fileb )

        return self.get_number_of_mutations(seq_a, seq_b)


    def main():

        fl1 = sys.argv[1]
        fl2 = sys.argv[2]

        a = get_seq_pdbfile(fl2)
        b = get_seq_pdbfile(fl1)

        seq_a = len(a)
        seq_b = len(b)

        print "Length of "+fl1+" is", seq_a
        print "Length of "+fl2+" is", seq_b
        print "Seq1 : "+fl2+" Seq2 : "+fl1
        tmp = print_def(a,b)

if __name__ == "__main__":
   analysemutations = AnalyseMutations()
   analysemutations.main()
