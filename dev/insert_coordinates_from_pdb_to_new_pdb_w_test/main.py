from numpy import mean,sqrt,var
import sys
from collections import defaultdict
from collections import OrderedDict
import os,shutil,sys,argparse

'''

The script is executed in the following way:

python ../main.py --file1 UM_10_M47H110H50C113_52H_-2_cui_21_10_2014_1_0001.pdb --file2 UM_48_M139H188H102C105_52H_-2_cui_21_10_2014_1_0001.pdb --residues 139,188,102,105 --include_ligand 1

and produces an output pdb, new.pdb. 

'''


class InsertCoordinates:

    def __init__(self):
        self.datafile1 = ""
        self.datafile2 = ""
        self.pdbone = []
        self.pdbtwo = []
        self.residues = []
        self.outputname = "new.pdb"
        self.previous_resnr = 0
        self.residuenumbers = ""
        self.include_ligand = 0
        self.hetatoms = []


    def generate_pdbfile(self):

        with open(self.datafile1,'r') as f:
            for line in f:
                self.pdbone.append( line )

        with open(self.datafile2,'r') as f:
            for line in f:
                if( line[23:26].strip() in self.residues ):
                    self.pdbtwo.append( line )
                elif( self.include_ligand != 0 and line[0:4] == "HETA" ):
                    self.hetatoms.append( line )

        with open(self.outputname,'w') as f:
            for line in self.pdbone:
                tmp_residue_nr = line[23:26].strip()
                if( tmp_residue_nr not in self.residues ):
                    f.write(line)

                elif( tmp_residue_nr != self.previous_resnr):
                    # print tmp_residue_nr, self.previous_resnr

                    for oline in self.pdbtwo:
                        if( oline[23:26].strip() == tmp_residue_nr ):
                            f.write(oline)
                    self.previous_resnr = tmp_residue_nr

            if(self.include_ligand != 0):
                for line in self.hetatoms:
                    f.write(line)


    def main(self):

        parser = argparse.ArgumentParser(description="Will insert coordinates from one pdb into another one")

        # get the initial rosetta design as input
        parser.add_argument("--file1", dest="datafile1", help=" PDB file 1" )

        parser.add_argument("--file2", dest="datafile2", help=" PDB file 2" )

        parser.add_argument("--residues", dest="residuenumbers", help="Residues that will be inserted into the pdbfile from datafile2" )

        parser.add_argument("--name", dest="outputname", help="Name of output file", default="new.pdb" )

        parser.add_argument("--include_ligand", dest="include_ligand", help="Include ligand", default=0 )


        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        # make list of residues
        self.residues = self.residuenumbers.split(",")

        self.generate_pdbfile()


if __name__ == "__main__":
   run = InsertCoordinates()
   run.main()


