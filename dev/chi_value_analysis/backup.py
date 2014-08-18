'''
ARG CHI1 N    CA   CB   CG
ARG CHI2 CA   CB   CG   CD
ARG CHI3 CB   CG   CD   NE
ARG CHI4 CG   CD   NE   CZ
ARG CHI5 CD   NE   CZ   NH1

ASN CHI1 N    CA   CB   CG
ASN CHI2 CA   CB   CG   OD1

ASP CHI1 N    CA   CB   CG
ASP CHI2 CA   CB   CG   OD1

CYS CHI1 N    CA   CB   SG

GLN CHI1 N    CA   CB   CG
GLN CHI2 CA   CB   CG   CD
GLN CHI3 CB   CG   CD   OE1

GLU CHI1 N    CA   CB   CG
GLU CHI2 CA   CB   CG   CD
GLU CHI3 CB   CG   CD   OE1

HIS CHI1 N    CA   CB   CG
HIS CHI2 CA   CB   CG   ND1

ILE CHI1 N    CA   CB   CG1
ILE CHI2 CA   CB   CG1  CD1

LEU CHI1N    CA   CB   CG
LEU CHI2CA   CB   CG   CD1

LYS CHI1 N    CA   CB   CG
LYS CHI2 CA   CB   CG   CD
LYS CHI3 CB   CG   CD   CE
LYS CHI4 CG   CD   CE   NZ

MET CHI1 N    CA   CB   CG
MET CHI2 CA   CB   CG   SD
MET CHI3 CB   CG   SD   CE

PHE CHI1 N    CA   CB   CG
PHE CHI2 CA   CB   CG   CD1

PRO CHI1 N    CA   CB   CG
PRO CHI2 CA   CB   CG   CD
PRO CHI3 CB   CG   CD   N
PRO CHI4 CG   CD   N    CA

SER CHI1N    CA   CB   OG

THR CHI1 N    CA   CB   OG1

TRP CHI1 N    CA   CB   CG
TRP CHI2 CA   CB   CG   CD1

TYR CHI1 N    CA   CB   CG
TYR CHI2 CA   CB   CG   CD1

VAL CHI1 N    CA   CB   CG1

( http://www.bchs.uh.edu/~glegge/phospho_lib/lib/tordef_phos.lib )


'''

#!/usr/bin/python
from optparse import OptionParser
from math import sqrt
from numpy import *
import os,shutil

'''

Chi-angle analysis. 

@Requires pdb file and list with residue numbers to be analysed e.g.,

python chi_analysis.py -f PDBFILE -l "1 2 3"

@parameter f: PDB file 
@type       : string

@parameter l: Three ligand names with Zn and its first ligand
@type       : string

'''

# Requires PDB files
# Returns list with lines in file
def get_pdbfile(pdbfile):
    fl = open(pdbfile,'r')
    pdb = fl.readlines()
    fl.close()
    return pdb

# Requires 4 vectors (numpy.array([]))
# Generates plane between four vectors
# Return angle between the two planes
def get_dihedral_angle(v1,v2,v3,v4):
    v12 = v2-v1
    v23 = v3-v2
    v34 = v4-v3
    # Normal vectors are generated
    n1 = cross(v12,v23)
    n2 = cross(v23,v34)

    angle = get_angle(n1,n2)
    # Getting the sign of the angle
    sign = dot(n1,v34)
    if sign < 0:
        angle = 360 - angle
    return angle

# Return angle between two vectors in degrees 
def get_angle(n1,n2):
    
    dot_p = dot(n1,n2) 
    
    norm_x = sqrt((n1*n1).sum())                                                                     
    norm_y = sqrt((n2*n2).sum()) 
    
    cos_angle = dot_p / norm_x / norm_y # cosinus of angle between x and y 
    
    angle = arccos(cos_angle)*57.3  
    
    return angle              


def get_residue_atoms(residue_name):
    residue_atom_names = {
        'ARG' : ['N','CA','CB','CG','CD','NE','CZ','NH1'],
        'ASN' : ['N','CA','CB','CG','OD1'],
        'ASP' : ['N','CA','CB','CG','OD1'],
        'CYS' : ['N','CA','CB','SG'],
        'GLN' : ['N','CA','CB','CG','CD','OE1'],
        'GLU' : ['N','CA','CB','CG','CD','OE1'],
        'HIS' : ['N','CA','CB','CG','ND1'],
        'ILE' : ['N','CA','CB','CG1','CD1'],
        'LEU' : ['N','CA','CB','CG','CD1'],
        'LYS' : ['N','CA','CB','CG','CD','CE','NZ'],
        'MET' : ['N','CA','CB','CG','SD','CE'],
        'PHE' : ['N','CA','CB','CG','CD1'],
        'PRO' : ['N','CA','CB','CG','CD'],
        'SER' : ['N','CA','CB','OG'],
        'THR' : ['N','CA','CB','OG1'],
        'TRP' : ['N','CA','CB','CG','CD1'],
        'TYR' : ['N','CA','CB','CG','CD1'],
        'VAL' : ['N','CA','CB','CG1']

        }


    return residue_atom_names[residue_name]


def get_chi_angles_residue(residue_number,pdbfile):
    residue_xyz_coordinates = []

    for line in pdbfile:

        if ( line[0:4] == 'ATOM' ):


            if ( line[22:26].strip() == residue_number ):

                residue_name = line[17:20]
                if ( residue_name in ['ALA', 'GLY'] ):
                    print "This residue type does not have chi angles: ", residue_name,residue_number
                    break
                    
                residue_atoms = get_residue_atoms(residue_name)
                
                if ( line[13:16].strip() in residue_atoms ) :

                    residue_xyz_coordinates.append(array(
                            [float(line[31:39]),
                             float(line[39:47]),
                             float(line[47:55])]))

    return residue_xyz_coordinates, residue_name       

def write_chi_angles(residue_xyz_coordinates):
    number_of_angles = len(residue_xyz_coordinates)
    torsion = []

    while number_of_angles >= 4:

        angle = get_dihedral_angle(residue_xyz_coordinates[0],
                                     residue_xyz_coordinates[1],
                                     residue_xyz_coordinates[2],
                                     residue_xyz_coordinates[3])
        if angle > 180:
            angle = angle - 360
        torsion.append(round(angle,2))
        residue_xyz_coordinates.pop(0)
        number_of_angles = number_of_angles - 1 

    return torsion


def write_chi_angles_pro(residue_xyz_coordinates):
    number_of_angles = len(residue_xyz_coordinates)

    torsion = []

    torsion.append(get_dihedral_angle(residue_xyz_coordinates[0],
                                residue_xyz_coordinates[1],
                                residue_xyz_coordinates[2],
                                residue_xyz_coordinates[3]))

    torsion.append(get_dihedral_angle(residue_xyz_coordinates[1],
                                residue_xyz_coordinates[2],
                                residue_xyz_coordinates[3],
                                residue_xyz_coordinates[4]))


    torsion.append(get_dihedral_angle(residue_xyz_coordinates[2],
                                residue_xyz_coordinates[3],
                                residue_xyz_coordinates[4],
                                residue_xyz_coordinates[0]))

    torsion.append(get_dihedral_angle(residue_xyz_coordinates[3],
                                residue_xyz_coordinates[4],
                                residue_xyz_coordinates[0],
                                residue_xyz_coordinates[1]))
    c_torsion = []

    for i in torsion:

        if i > 180:
            i = i - 360

        c_torsion.append(round(i,2))
    
    return c_torsion



def main():

    parser = OptionParser()

    parser.add_option('-f',dest='PDB',
                      help='Cleaned pdb file with metal ion present default=ZN')

    parser.add_option('-l',dest='residue_number',default='25',
                      help='Residues to analyse for chi-values')


    (options,args) = parser.parse_args()
    
    pdbfile = get_pdbfile(options.PDB)

    # Open file to write to 

    # Loop over inout here
    
    residue_number = options.residue_number.split()
    
    for i in residue_number:
    
        vectors,residue_name = get_chi_angles_residue(i,pdbfile)

        if residue_name != 'PRO':

            chi_angles = write_chi_angles(vectors)

            print residue_name, chi_angles

        else:

            chi_angles = write_chi_angles_pro(vectors)

            print residue_name, chi_angles



if __name__ == "__main__":
    main()
