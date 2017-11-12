#!/usr/bin/env python
import sys
from numpy import *

"""

"""


# @requires pdbfile
# @return list with coordinates
def get_pdbfile(pdbfile):
    tmp_file = open(pdbfile,'r')
    pdbfile = tmp_file.readlines()
    tmp_file.close()
    return pdbfile


# replace residue name 
def replace_residue_name(new_residue_name,old_residue_name,pdbfile):
    new_pdbfile = []
    for line in pdbfile:
        if line[17:20] == old_residue_name:
            line = line[0:17]+new_residue_name+line[20:]
        new_pdbfile.append(line)
    return new_pdbfile
    

def write_file(new_pdbfile):
    new_file = open('correct_a.pdb','w')
    for line in new_pdbfile:
        new_file.write(line)


def replace_hetatom_names(old_pdbfile):
    phosphate_part = 'P2M'
    base = 'X5 '
    new_pdbfile = []

    protein = []
    phosphate = ['']*5
    base_coor = ['']*18


    for line in old_pdbfile:

        if line[0:4] == 'HETA':
            
            if line[13:16].strip() == 'P1':
                line = line[0:13]+'P  '+line[16:17]+phosphate_part+line[20:]
                phosphate[1] = line

            elif line[13:16].strip() == 'C1':
                line = line[0:13]+'C5\''+line[16:17]+base+line[20:]
                base_coor[0] = line
                
            elif line[13:16].strip() == 'N1':
                line = line[0:13]+'N3 '+line[16:17]+base+line[20:]
                base_coor[15] = line

            elif line[13:16].strip() == 'O7':
                line = line[0:13]+'OP3'+line[16:17]+phosphate_part+line[20:]
                phosphate[4] = line 

            elif line[13:16].strip() == 'C2':
                line = line[0:13]+'C4\''+line[16:17]+base+line[20:]
                base_coor[1] = line 

            elif line[13:16].strip() == 'N2':
                line = line[0:13]+'N4 '+line[16:17]+base+line[20:]
                base_coor[14] = line 

            elif line[13:16].strip() == 'O5':
                line = line[0:13]+'OP2'+line[16:17]+phosphate_part+line[20:]
                phosphate[3] = line 

            elif line[13:16].strip() == 'C3':
                line = line[0:13]+'C3\''+line[16:17]+base+line[20:]
                base_coor[2] = line 

            elif line[13:16].strip() == 'N3':
                line = line[0:13]+'N5 '+line[16:17]+base+line[20:]
                base_coor[10] = line 

            elif line[13:16].strip() == 'O4':
                line = line[0:13]+'OP1'+line[16:17]+phosphate_part+line[20:]
                phosphate[2] = line 

            elif line[13:16].strip() == 'C4':
                line = line[0:13]+'C2\''+line[16:17]+base+line[20:]
                base_coor[4] = line 

            elif line[13:16].strip() == 'O1':
                line = line[0:13]+'OR '+line[16:17]+phosphate_part+line[20:]
                phosphate[0] = line 

            elif line[13:16].strip() == 'C5':
                line = line[0:13]+'C1\''+line[16:17]+base+line[20:]
                base_coor[5] = line 

            elif line[13:16].strip() == 'O2':
                line = line[0:13]+'O4\''+line[16:17]+base+line[20:]
                base_coor[6] = line 

            elif line[13:16].strip() == 'C6':
                line = line[0:13]+'C4 '+line[16:17]+base+line[20:]
                base_coor[13] = line 

            elif line[13:16].strip() == 'O3':
                line = line[0:13]+'O3\''+line[16:17]+base+line[20:]
                base_coor[3] = line 

            elif line[13:16].strip() == 'C7':
                line = line[0:13]+'C5 '+line[16:17]+base+line[20:]
                base_coor[9] = line 

            elif line[13:16].strip() == 'O6':
                line = line[0:13]+'O2 '+line[16:17]+base+line[20:]
                base_coor[17] = line 

            elif line[13:16].strip() == 'C8':
                line = line[0:13]+'C2 '+line[16:17]+base+line[20:]
                base_coor[16] = line 

            elif line[13:16].strip() == 'O8':
                line = line[0:13]+'O51'+line[16:17]+base+line[20:]
                base_coor[11] = line 

            elif line[13:16].strip() == 'C9':
                line = line[0:13]+'C6 '+line[16:17]+base+line[20:]
                base_coor[8] = line 

            elif line[13:16].strip() == 'O9':
                line = line[0:13]+'O52'+line[16:17]+base+line[20:]
                base_coor[12] = line 

            elif line[13:16].strip() == 'C10':
                line = line[0:13]+'C1 '+line[16:17]+base+line[20:]
                base_coor[7] = line 
        elif line[0:4] == 'ATOM':
            protein.append(line)

    new_pdbfile = protein+phosphate+base_coor
    return new_pdbfile

        
def main():
    old_pdbfile = get_pdbfile('a.pdb')
    new_pdbfile = replace_residue_name('DZM','ZMP',old_pdbfile)

    new_pdbfile = replace_hetatom_names(new_pdbfile)

    write_file(new_pdbfile)
        
if __name__ == "__main__":
    main()
