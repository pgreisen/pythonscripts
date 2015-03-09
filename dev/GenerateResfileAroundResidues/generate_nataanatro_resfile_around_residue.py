#!/usr/bin/env python
import sys
from numpy import *
import argparse

"""
Generate resfile for enzyme design

@requires pdb file with ligand present
@return resfile with natro for all residue > 10 AA and nataa > 8 AA

"""


# return coordinates of residue
def get_ligand_coordinates(pdbfile,residuenumber):
    pdbfl = open(pdbfile,'r')
    ligand_coor = []
    for line in pdbfl:
        if line[22:26] == residuenumber:
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
    # print start,end
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
    DISTANCE2 = 10

    # File name of pdb file
    pdbfile = sys.argv[1]

    residuenumber = sys.argv[2]

    # Collecting the ligand coordinates in
    # a vector

    ligand = get_ligand_coordinates(pdbfile,residuenumber)
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
                    
                
                if tmp_length > DISTANCE1 and tmp_length < DISTANCE2:
                    
                    pdb_resid = int(str(line[23:26]).rstrip())
                    nataa.append(pdb_resid)


                elif tmp_length > DISTANCE2:

                    pdb_resid = int(str(line[23:26]).rstrip())
                    natro.append(pdb_resid)
                    
                    
                elif tmp_length < DISTANCE1:

                    pdb_resid = int(str(line[23:26]).rstrip())
                    design.append(pdb_resid)

                else:
                    print 'Bug in program'


    # get the rigth amino acids
    nataa = set(nataa).difference(set(design))
    natro = set(natro).difference(set(nataa))
    natro = set(natro).difference(set(design)) 

    # Get constraint residues
    cst_residues = get_constraint_residues(pdbfl)

    start,end = get_length_protein(pdbfile)

    # resfile
    rs_file = open('resfile','w')
    rs_file.write('start\n')
    
    rs_1 = ' A NATAA\n'
    rs_2 = ' A NATRO\n'
    
    while start <= end:
        
        if start in natro:
            rs_file.write('\t'+str(start)+rs_2)

        elif start in cst_residues:
            rs_file.write('\t'+str(start)+rs_2)

        else:
            rs_file.write('\t'+str(start)+rs_1)


            
        start = start + 1
        
if __name__ == "__main__":
    main()
