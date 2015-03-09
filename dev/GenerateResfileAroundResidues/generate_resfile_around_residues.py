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
        self.DISTANCE1 = 6
        self.DISTANCE2 = 8
        self.residuenumber = ""
        self.residues = []
        self.pdbfile = []
        self.datafile = []


    def get_pdbfile(self,file):
        with open(file,'r') as f:
            for line in f:
                self.pdbfile.append( line )


    # return coordinates of residue
    # get_residue_coordinates
    def get_residue_coordinates(self):
        # Debug

        ligand_coor = []
        for line in self.pdbfile:
            if line[22:26].strip() in self.residues:
                x = str(line[30:38]).rstrip()
                y = str(line[38:46]).rstrip()
                z = str(line[46:54]).rstrip()
                ligand_coor.append(array([float(x),float(y),float(z)]))
        print "Ligand Coordinates ", ligand_coor
        return ligand_coor

    def get_length_protein(self):
        first = '0'
        start = ''
        end = ''
        for line in self.pdbfile:
            if line[0:4] == 'ATOM':
                if first == '0':
                    start = str(line[23:26]).rstrip()
                    first = '1'
                elif first == '1':
                    end = str(line[23:26]).rstrip()
        return int(start),int(end)

    # @requires pdbfile
    # @return list with coordinates
    def get_pdbfile(self, pdbfile):
        tmp_file = open(pdbfile,'r')
        for line in tmp_file:
            self.pdbfile.append(line)
        tmp_file.close()
            
    # @requires pdblist
    # @return list with constraint residues
    def get_constraint_residues(pdblist):
        cst = []
        for line in pdblist:
            if line[0:4] == 'REMA':
                cst.append(int(line.split()[11]))
            
            else:
                break
        return set(cst)

    def main(self):
        parser = argparse.ArgumentParser(description="Generate Residue File Around A List of Residues")
        parser.add_argument('-f',dest='datafile', help='PDB files' )
        parser.add_argument('--residuenumber',dest='residuenumber', help='Residues with PDB numbering which are used to generate the resfile' )
        parser.add_argument('--distance1',dest='distance1', default=6,type=int, help="If atoms are below this value the residue is allowed to repack")
        parser.add_argument('--distance2',dest='distance2', default=8,type=int, help="Above this threshold residues are not allowed to repack")
        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        print self.residuenumber
        self.residues = self.residuenumber.split(",")

        # Get pdbfile
        self.get_pdbfile( self.datafile)


        ligand = self.get_residue_coordinates()
    
        # Collecting residue id for residues within or outside
        # ligand
        nataa = []
        natro = []
        # We do not added any design positions
        # will be done later
        # design = []

        for i in ligand:

            for line in self.pdbfile:

                if line[0:4] == 'ATOM':
                    x = str(line[30:38]).rstrip()
                    y = str(line[38:46]).rstrip()
                    z = str(line[46:54]).rstrip()

                    tmp_vector = array([float(x),float(y),float(z)])

                    tmp_length = linalg.norm( tmp_vector - i)

                    pdb_resid = int(str(line[23:26]).rstrip())
                
                    if tmp_length > self.DISTANCE1:
                        pdb_resid = int(str(line[23:26]).rstrip())
                        natro.append(pdb_resid)

                    elif tmp_length <= self.DISTANCE1:
                        pdb_resid = int(str(line[23:26]).rstrip())
                        nataa.append(pdb_resid)

                    else:
                        print 'Bug in program'

        #
        # print nataa
        nataa = set(nataa) #.difference(set(natro))

        natro = set(natro).difference(set(nataa))

        # TODO add cst option to script
        # Get constraint residues
        # cst_residues = get_constraint_residues(pdbfl)

        start,end = self.get_length_protein()

        # resfile
        rs_file = open('resfile','w')
        rs_file.write('start\n')
    
        rs_1 = ' A NATAA\n'
        rs_2 = ' A NATRO\n'
    
        while start <= end:
        
            if start in natro:
                rs_file.write('\t'+str(start)+rs_2)
            elif( start in nataa ):
                rs_file.write('\t'+str(start)+rs_1)

            start = start + 1
        
if __name__ == "__main__":
    run = GenerateResfileAroundResidues()
    run.main()
