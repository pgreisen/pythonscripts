#!/usr/bin/env python
import sys
from numpy import *
import argparse

"""
Generate resfile for enzyme design

@requires pdb file with ligand present
@return resfile with natro for all residue > 10 AA and nataa > 8 AA

"""

class GenerateResfileAroundResidues:


    def __init__(self):
        self.pdbfile = ""
        self.backboneatoms = ['N','CA','C','O']


    def write_poly_glycine_pdbfile(self,file):

        polyglycinepdb = open("polyglycine.pdb",'w')


        with open(file,'r') as f:
            for line in f:

                tmp = line.split()

                if( len(tmp) > 7 and tmp[0] == "ATOM" and tmp[2] in self.backboneatoms ):
                    polyglycinepdb.write(line[0:17]+"GLY"+line[20:])

        polyglycinepdb.close()


    def main(self):

        parser = argparse.ArgumentParser(description="Generate Residue File Around A List of Residues")
        parser.add_argument('-f',dest='pdbfile', help='PDB files' )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])


        self.write_poly_glycine_pdbfile(self.pdbfile)

        
if __name__ == "__main__":
    run = GenerateResfileAroundResidues()
    run.main()
