#!/usr/bin/env python
'''
Computes the rmsd difference between two poses per residue

python get_rmsd_per_residue --file1 PDB1 --file2 PDB2


RMSD = sum(N) = sqrt([(xi - xj)**2] / N )

'''

import argparse
from pdbfile import *

class GetRmsdPerResidue:


    def __init__(self):
        self.rmsd_per_residue = {}
        self.counts = {}


    def get_rmsd_between_residues(self, p1_coordinates, p2_coordinates):
        rmsd = {}
        tmp = 0
        for key in p1_coordinates:
            new_key = key[10:19]
            if(key[10:19] in self.rmsd_per_residue):
                self.rmsd_per_residue[new_key] += linalg.norm((p1_coordinates[key] - p2_coordinates[key]))
                self.counts[new_key] += 1

        for key in self.rmsd_per_residue:
            rmsd[key] = sqrt(self.rmsd_per_residue[key] / self.counts[key])

        return rmsd



    def get_coordinates(self,pdbfile):
        pdbcoordinates = {}
        rmsd_per_residue = {}
        for line in pdbfile:
            if(line[0:4] == "ATOM" and line[13:14] != 'H'):
                pdbcoordinates[str(line[7:26])] = array([float(line[30:39]), float(line[39:47]), float(line[47:55])])
                key2 = str(line[17:26])
                self.rmsd_per_residue[key2] = 0
                self.counts[key2] = 0
        return pdbcoordinates




    def main(self):

        parser = argparse.ArgumentParser(description=" Computes the RMSD between two poses and plots it as well as dumping a rmsd.dat ")
        parser.add_argument("--file1", dest="pdbfile1", help="First pdb file")
        parser.add_argument("--file2", dest="pdbfile2", help="First pdb file")

        input_variables = parser.parse_args()
        # initialize pdbfile
        pf = pdbfile()

        pdbfile1 = pf.read_file(input_variables.pdbfile1)
        pdbfile2 = pf.read_file(input_variables.pdbfile2)

        # Coordinates are stored in dictionary
        p1_coordinates = self.get_coordinates(pdbfile1)
        p2_coordinates = self.get_coordinates(pdbfile2)

        rmsd = self.get_rmsd_between_residues(p1_coordinates,p2_coordinates)

        rmsd_file = open("rmsd.dat","w")

        for key in rmsd:
            rmsd_file.write(key+"\t"+str(round(rmsd[key],3))+"\n")
        rmsd_file.close()

    def test(self):

        pdbfile1 = pf.read("2RIN_chainA_0001.pdb")
        pdbfile2 = pf.read("docking_2RIN_chainA_0001_0026.pdb")


if __name__ == "__main__":
   run = GetRmsdPerResidue()
   run.main()