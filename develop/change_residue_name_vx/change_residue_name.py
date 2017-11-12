#!/usr/bin/env python
import sys
from numpy import *

"""
Instead of all the elif make a dictionary as input


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
    

def write_file(new_pdbfile,pdbname):
    new_file = open('c_'+pdbname,'w')
    for line in new_pdbfile:
        new_file.write(line)

# Takes to lists of atoms
def replace_hetatom_names(correct_list,replace_list):
    
    # Collect the substituted atom and residue names
    # Need to add OH- to replace atoms
    # Need to add Zn atoms etc later on
    replace_atomname_residuename = ['']*50
    length_of_list = len(replace_list)
    for i in range(length_of_list):
        # P1 -> P1

        if(replace_list[i][13:16].strip() == 'P1'):

            line = correct_list[20][0:26]+replace_list[i][26:55]+correct_list[20][55:]

            replace_atomname_residuename[20] = line

        # O2 -> O4 
        elif(replace_list[i][13:16].strip() == 'O2'):

            line = correct_list[4][0:26]+replace_list[i][26:55]+correct_list[4][55:]

            replace_atomname_residuename[4] = line

        # S1 -> S1
        elif(replace_list[i][13:16].strip() == 'S1'):

            line = correct_list[21][0:26]+replace_list[i][26:55]+correct_list[21][55:]

            replace_atomname_residuename[21] = line

        # O1 -> O3
        elif(replace_list[i][13:16].strip() == 'O1'):

            line = correct_list[17][0:26]+replace_list[i][26:55]+correct_list[17][55:]

            replace_atomname_residuename[17] = line
        # N1 -> N1
        elif(replace_list[i][13:16].strip() == 'N1'):

            line = correct_list[5][0:26]+replace_list[i][26:55]+correct_list[5][55:]

            replace_atomname_residuename[5] = line
        # C9 -> C10
        elif(replace_list[i][13:16].strip() == 'C9'):

            line = correct_list[14][0:26]+replace_list[i][26:55]+correct_list[14][55:]

            replace_atomname_residuename[14] = line
        # C10 -> C11
        elif(replace_list[i][13:16].strip() == 'C10'):

            line = correct_list[15][0:26]+replace_list[i][26:55]+correct_list[15][55:]

            replace_atomname_residuename[15] = line

        # C11 -> C12
        elif(replace_list[i][13:16].strip() == 'C11'):

            line = correct_list[16][0:26]+replace_list[i][26:55]+correct_list[16][55:]

            replace_atomname_residuename[16] = line

        # C4 -> C5
        elif(replace_list[i][13:16].strip() == 'C4'):

            line = correct_list[9][0:26]+replace_list[i][26:55]+correct_list[9][55:]

            replace_atomname_residuename[9] = line
        # C3 -> C4
        elif(replace_list[i][13:16].strip() == 'C3'):

            line = correct_list[8][0:26]+replace_list[i][26:55]+correct_list[8][55:]

            replace_atomname_residuename[8] = line
        # C2 - C3
        elif(replace_list[i][13:16].strip() == 'C2'):

            line = correct_list[7][0:26]+replace_list[i][26:55]+correct_list[7][55:]

            replace_atomname_residuename[7] = line
        # C8 -> C9
        elif(replace_list[i][13:16].strip() == 'C8'):

            line = correct_list[13][0:26]+replace_list[i][26:55]+correct_list[13][55:]

            replace_atomname_residuename[13] = line
        # C7 -> C8
        elif(replace_list[i][13:16].strip() == 'C7'):

            line = correct_list[12][0:26]+replace_list[i][26:55]+correct_list[12][55:]

            replace_atomname_residuename[12] = line

        # C1 -> C2
        elif(replace_list[i][13:16].strip() == 'C1'):

            line = correct_list[6][0:26]+replace_list[i][26:55]+correct_list[6][55:]

            replace_atomname_residuename[6] = line

        # C5 -> C6
        elif(replace_list[i][13:16].strip() == 'C5'):

            line = correct_list[10][0:26]+replace_list[i][26:55]+correct_list[10][55:]

            replace_atomname_residuename[10] = line

        # C6 -> C7
        elif(replace_list[i][13:16].strip() == 'C6'):

            line = correct_list[11][0:26]+replace_list[i][26:55]+correct_list[11][55:]

            replace_atomname_residuename[11] = line

        # H1
        elif(replace_list[i][13:16].strip() == 'H1'):

            line = correct_list[22][0:26]+replace_list[i][26:55]+correct_list[22][55:]

            replace_atomname_residuename[22] = line
        # H2
        elif(replace_list[i][13:16].strip() == 'H2'):

            line = correct_list[23][0:26]+replace_list[i][26:55]+correct_list[23][55:]

            replace_atomname_residuename[23] = line
        # H3
        elif(replace_list[i][13:16].strip() == 'H3'):

            line = correct_list[24][0:26]+replace_list[i][26:55]+correct_list[24][55:]

            replace_atomname_residuename[24] = line

        # H4
        elif(replace_list[i][13:16].strip() == 'H4'):

            line = correct_list[25][0:26]+replace_list[i][26:55]+correct_list[25][55:]

            replace_atomname_residuename[25] = line


        # H5 
        elif(replace_list[i][13:16].strip() == 'H5'):

            line = correct_list[26][0:26]+replace_list[i][26:55]+correct_list[26][55:]

            replace_atomname_residuename[26] = line
        # H6
        elif(replace_list[i][13:16].strip() == 'H6'):

            line = correct_list[27][0:26]+replace_list[i][26:55]+correct_list[27][55:]

            replace_atomname_residuename[27] = line

        # H7
        elif(replace_list[i][13:16].strip() == 'H7'):

            line = correct_list[28][0:26]+replace_list[i][26:55]+correct_list[28][55:]

            replace_atomname_residuename[28] = line
        # H8
        elif(replace_list[i][13:16].strip() == 'H8'):

            line = correct_list[29][0:26]+replace_list[i][26:55]+correct_list[29][55:]

            replace_atomname_residuename[29] = line
        # H9
        elif(replace_list[i][13:16].strip() == 'H9'):

            line = correct_list[30][0:26]+replace_list[i][26:55]+correct_list[30][55:]

            replace_atomname_residuename[30] = line

        # H10
        elif(replace_list[i][13:16].strip() == 'H10'):

            line = correct_list[31][0:26]+replace_list[i][26:55]+correct_list[31][55:]

            replace_atomname_residuename[31] = line
        # H11
        elif(replace_list[i][13:16].strip() == 'H11'):

            line = correct_list[32][0:26]+replace_list[i][26:55]+correct_list[32][55:]

            replace_atomname_residuename[32] = line
        # H12
        elif(replace_list[i][13:16].strip() == 'H12'):

            line = correct_list[33][0:26]+replace_list[i][26:55]+correct_list[33][55:]

            replace_atomname_residuename[33] = line

        # H13
        elif(replace_list[i][13:16].strip() == 'H13'):

            line = correct_list[34][0:26]+replace_list[i][26:55]+correct_list[34][55:]

            replace_atomname_residuename[34] = line
        # H14
        elif(replace_list[i][13:16].strip() == 'H14'):

            line = correct_list[35][0:26]+replace_list[i][26:55]+correct_list[35][55:]

            replace_atomname_residuename[35] = line
        # H15
        elif(replace_list[i][13:16].strip() == 'H15'):

            line = correct_list[36][0:26]+replace_list[i][26:55]+correct_list[36][55:]

            replace_atomname_residuename[36] = line

        # H16
        elif(replace_list[i][13:16].strip() == 'H16'):

            line = correct_list[37][0:26]+replace_list[i][26:55]+correct_list[37][55:]

            replace_atomname_residuename[37] = line
        # H17
        elif(replace_list[i][13:16].strip() == 'H17'):

            line = correct_list[38][0:26]+replace_list[i][26:55]+correct_list[38][55:]

            replace_atomname_residuename[38] = line

        # H18
        elif(replace_list[i][13:16].strip() == 'H18'):

            line = correct_list[39][0:26]+replace_list[i][26:55]+correct_list[39][55:]

            replace_atomname_residuename[39] = line

        # H19
        elif(replace_list[i][13:16].strip() == 'H19'):

            line = correct_list[40][0:26]+replace_list[i][26:55]+correct_list[40][55:]

            replace_atomname_residuename[40] = line

        # H20
        elif(replace_list[i][13:16].strip() == 'H20'):

            line = correct_list[41][0:26]+replace_list[i][26:55]+correct_list[41][55:]

            replace_atomname_residuename[41] = line

        # H21
        elif(replace_list[i][13:16].strip() == 'H21'):

            line = correct_list[42][0:26]+replace_list[i][26:55]+correct_list[42][55:]

            replace_atomname_residuename[42] = line


        # H22
        elif(replace_list[i][13:16].strip() == 'H22'):

            line = correct_list[43][0:26]+replace_list[i][26:55]+correct_list[43][55:]

            replace_atomname_residuename[43] = line


        elif(replace_list[i][13:16].strip() == 'H23'):

            line = correct_list[44][0:26]+replace_list[i][26:55]+correct_list[44][55:]

            replace_atomname_residuename[44] = line
        # H24
        elif(replace_list[i][13:16].strip() == 'H24'):

            line = correct_list[45][0:26]+replace_list[i][26:55]+correct_list[45][55:]

            replace_atomname_residuename[45] = line
            
        # H25
        elif(replace_list[i][13:16].strip() == 'H25'):

            line = correct_list[46][0:26]+replace_list[i][26:55]+correct_list[46][55:]

            replace_atomname_residuename[46] = line
        # H26
        elif(replace_list[i][13:16].strip() == 'H26'):

            line = correct_list[47][0:26]+replace_list[i][26:55]+correct_list[47][55:]

            replace_atomname_residuename[47] = line
        # H27
        elif(replace_list[i][13:16].strip() == 'H27'):

            line = correct_list[48][0:26]+replace_list[i][26:55]+correct_list[48][55:]

            replace_atomname_residuename[48] = line

    replace_atomname_residuename[0] = correct_list[0]
    replace_atomname_residuename[1] = correct_list[1]
    replace_atomname_residuename[2] = correct_list[2]
    replace_atomname_residuename[3] = correct_list[3]
    replace_atomname_residuename[4] = correct_list[4]
    replace_atomname_residuename[18] = correct_list[18]
    replace_atomname_residuename[19] = correct_list[19]
    replace_atomname_residuename[49] = correct_list[49]

    return replace_atomname_residuename

        
def main():
    original = sys.argv[1]
    model = sys.argv[2]

    template_pdbfile = get_pdbfile(original)
    model_pdbfile = get_pdbfile(model)

    new_pdbfile = replace_hetatom_names(template_pdbfile,model_pdbfile)

    write_file(new_pdbfile,model)
        
if __name__ == "__main__":
    main()
