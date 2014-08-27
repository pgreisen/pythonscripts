#!/usr/bin/env python
import sys
from numpy import *

"""
Generate resfile for design

@requires pdb file with ligand present
@return resfile with natro for all residue > 8 AA and nataa > 6 AA and < 8 AA.

"""

def get_bigger_amino_acids(design,pdbfile):
    design_shell_residues = {
        'ALA' : 'AVLIFM',
        'VAL' : 'VILMF',
        'SER' : 'SNQT',
        'ASN' : 'NQ',
        'ILE' : 'ILFMYW',
        'LEU' : 'ILFMYW'
        }


    bigger_amino_acids = {}
    for line in pdbfile:

        if line[0:4] == 'ATOM':

            residuenumber = int(str(line[23:26]).rstrip())

            if residuenumber in design:

                residuename = line[17:20]

                for key, value in design_shell_residues.iteritems():

                    if key == residuename:
                        bigger_amino_acids[residuenumber] = value

    return bigger_amino_acids
    

    


# @require pdbfile and atom specifier
# @return list with protein or ligand

def get_ligand_coordinates(pdbfile,chemical='HETA'):
    pdbfl = open(pdbfile,'r')
    ligand_coor = []
    for line in pdbfl:
        if line[0:4] == chemical:
            x = str(line[30:38]).rstrip()
            y = str(line[38:46]).rstrip()
            z = str(line[46:54]).rstrip()
            ligand_coor.append(array([float(x),float(y),float(z)]))
    return ligand_coor

def get_length_protein(pdbfile):
    pdbfl = open(pdbfile,'r')
    first = '0'
    start = ''
    end = ''
    for line in pdbfl:
        if line[0:4] == 'ATOM':
            if first == '0':
                start = str(line[23:26]).rstrip()
                first = '1'
            elif first == '1':
                end = str(line[23:26]).rstrip()
    return int(start),int(end)
            

# @requires pdbfile
# @return list with coordinates

def get_pdbfile(pdbfile):
    pdbfl = []
    tmp_file = open(pdbfile,'r')
    for line in tmp_file:
        pdbfl.append(line)
    tmp_file.close()
    return pdbfl
            
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


def main():

    # Distance which the residue should be left alone
    DISTANCE1 = 6
    DISTANCE2 = 8

    # File name of pdb file
    pdbfile = sys.argv[1]

    # Collecting the ligand coordinates in
    # a vector
    ligand = get_ligand_coordinates(pdbfile)
    # print 'Length of ligand', len(ligand)
    
    # Collecting residue id for residues within or outside
    # ligand
    nataa = []
    natro = []
    design = []

    # Get list of pdb-file
    pdbfl = get_pdbfile(pdbfile)

    
    for i in ligand:

        for line in pdbfl:

            if line[0:4] == 'ATOM':
                
                x = str(line[30:38]).rstrip()
                y = str(line[38:46]).rstrip()
                z = str(line[46:54]).rstrip()

                tmp_vector = array([float(x),float(y),float(z)])

                tmp_length = linalg.norm(tmp_vector - i)

                pdb_resid = int(str(line[23:26]).rstrip())
                    
                
                if tmp_length > DISTANCE1 and tmp_length <= DISTANCE2:
                    
                    pdb_resid = int(str(line[23:26]).rstrip())
                    nataa.append(pdb_resid)


                elif tmp_length > DISTANCE2:
                    
                    pdb_resid = int(str(line[23:26]).rstrip())
                    natro.append(pdb_resid)
                    print pdb_resid
                    
                elif tmp_length <= DISTANCE1:

                    pdb_resid = int(str(line[23:26]).rstrip())
                    design.append(pdb_resid)

                else:
                    print 'Bug in program'


    # get the rigth amino acids
    # Designable residue
    design = set(design)
    # Get the difference between design and nataa which should correspond 
    # to amino acids with a distance larger than 6 AA
    nataa = set(nataa).difference(set(design))
    # Get the difference between nataa and natro which should correspond 
    # to amino acids with a distance larger than 8 AA
    natro = set(natro).difference(set(nataa))
    # Remember some of the designed residues can have atoms with 
    # distance longer than 8 AA
    natro = set(natro).difference(set(design))
    # Get constraint residues
    cst_residues = get_constraint_residues(pdbfl)

    start,end = get_length_protein(pdbfile)

    # resfile
    rs_file = open('resfile','w')
    rs_file.write('start\n')
    
    rs_1 = ' A NATAA\n'
    rs_2 = ' A NATRO\n'

    rs_3 = ' A PIKAA '

    extend_amino_acids = get_bigger_amino_acids(design,pdbfl)
    print natro
    # Debug 
    # for k,v in extend_amino_acids.iteritems():
    #    print k,v
    #
    print design
    while start <= end:
        
        if start in natro:
            rs_file.write('\t'+str(start)+rs_2)

        elif start in cst_residues:
            rs_file.write('\t'+str(start)+rs_2)

        if start in design:

            if int(start) in extend_amino_acids.keys():

                rs_file.write('\t'+str(start)+rs_3+extend_amino_acids[start]+'\n')
            else:
                rs_file.write('\t'+str(start)+rs_1)

        else:
            rs_file.write('\t'+str(start)+rs_1)

        start = start + 1


        
if __name__ == "__main__":
    main()
