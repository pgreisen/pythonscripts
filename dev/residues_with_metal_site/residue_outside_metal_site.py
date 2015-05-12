#!/usr/bin/python
from optparse import OptionParser
from math import sqrt
from numpy import *
import os,shutil

'''
Fix_pdbfile get chain A as well as hetatms connected to chain A.
It takes all a pdb file as input and returns
the same file.

@parameter f: PDB file with metal site
@type       : string
@parameter m: Metal ion - default ZN 
@type       : string
@parameter n: Name of rosetta output pdb file 
@type       : string
@parameter l: Three ligand names with Zn and its first ligand
@type       : string
@parameter a: Coordinate file of aligned ligand and is append to the Rosetta input
@type       : string
@return     : constraint file, pdb file to use as input for rosetta


The default metal-protein ligand distance is set to 2.9 AA and for non-protein it is set to 2.7 AA


'''


class MetalSiteGeometry:

    def __init__(self):
        self.distance = 10
        self.residues = []
        self.METAL = "ZN"
        self.metal_coor = {}
        self.ddg_file = "ddG_rosetta.dat"

        self.aa = {
            'H' :'HIS',
            'E' :'GLU',
            'D' :'ASP',
            'C' :'CYS',
            'K' :'LYS',
            'T' :'THR',
            'S' :'SER',
            'R' :'ARG',
            'Q' :'GLN',
            'N' :'ASN',
            'Y' :'TYR',
            'M' :'MET',
            'P' :'PRO',
            'A' :'ALA',
            'V' :'VAL',
            'L' :'LEU',
            'I' :'ILE',
            'F' :'PHE',
            'W' :'TRP',
            'G' :'GLY'
        }



    # Requires PDB files
    # Returns list with lines in file
    # and sets the metal ions for later geometry determination
    def get_pdbfile(self,pdbfile):
        pdb = []
        with open(pdbfile) as f:
            for line in f:
                if(line[0:4] == "ATOM"):
                    pdb.append( line )
                elif(line[0:4] == "HETA"):
                    atom_name = str(line[12:15]).strip()
                    het_id = line[13:14]
                    if atom_name == self.METAL:
                        # key is the residue name with the chain and residue number information
                        self.metal_coor[line[18:26]] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])

                    pdb.append(line)
        return pdb

    # Requires list
    # Returns one list
    def merge(self,seq):
        merged = []
        for s in seq:
            for x in s:
                merged.append(x)
        return merged


    # Requires pdb file
    # Returns protein atoms
    def get_atoms_pdb(self,PDB):
        atm = []
        for line in PDB:
            if line[0:4] == 'ATOM':
                atm.append(line)
        return atm

    # Requires pdb file
    # Returns heteroatoms
    def get_heteroatoms_pdb(self, PDB):
        atm = []
        for line in PDB:
            if line[0:4] == 'HETA':
                atm.append(line)
        return atm


    # Requires list with pdb lines, metal coordinates ( set in the constructor )
    # Returns dictionary with ligands < DISTANCE from metal ion
    def get_protein_ligand_metal(self, PDB ):

        # loop over metal sites
        for metal_site in self.metal_coor:

            metal_vec = self.metal_coor[metal_site]

            for line in PDB:

                res = line[17:20].rstrip()

                vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                # distance between the metal site and protein atom
                ds_lig = linalg.norm(vec-metal_vec)

                if ds_lig < self.distance:

                    residue_number = line[22:26].strip()

                    self.residues.append( residue_number )



    # Requires PDB
    # Returns chain of first pdb
    def get_chain(self,pdb):
        for line in pdb:
            if line[0:4] == 'ATOM':
                chain = line[21:22]
                break
        return str(chain)


    def get_list_of_substitutions(self):

        tmpfile = open("make_substitutions.sh",'w')

        with open(self.ddg_file,'r') as f:
            for line in f:
                tmp = line.split(',')

                native, residuenr, substitution = tmp[0].split('_')

                if residuenr not in self.residues:
                    #if ( substitution == 'P'):
                    #
                    #    print line
                    # print line, self.aa[substitution]

                    tmpfile.write("python ~/pythonscripts/dev/replace_amino_acids_at_position/replace_aa_at_position.py --aa "+self.aa[substitution]+" --position "+str(residuenr)+" -f $1\n")
        tmpfile.close()




    def main(self):
        parser = OptionParser()

        parser.add_option('-f',dest='PDB', help='PDB file with metal ion present default=ZN')
        parser.add_option('-d',dest='distance', help='Cutoff distance from the metal site')

        parser.add_option('--ddg',dest='ddg_file', help='File with results from the ddG application in Rosetta')


        (options, args) = parser.parse_args()

        PDB = self.get_pdbfile(options.PDB)

        self.get_protein_ligand_metal( PDB )

        self.residues = set( self.residues )

        self.get_list_of_substitutions()



if __name__ == "__main__":
    run = MetalSiteGeometry()
    run.main()
