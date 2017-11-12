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
            l_nr = 1 -1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # O2 -> O2
        elif(replace_list[i][13:16].strip() == 'O2'):

            l_nr = 7-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # S1 -> S1
        elif(replace_list[i][13:16].strip() == 'S1'):


            l_nr = 5-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # O1 -> O1
        elif(replace_list[i][13:16].strip() == 'O1'):

            l_nr = 4-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # N1 -> N1
        elif(replace_list[i][13:16].strip() == 'N1'):

            l_nr = 3-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # C9 -> C9
        elif(replace_list[i][13:16].strip() == 'C9'):
            l_nr = 14-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # C10 -> C10
        elif(replace_list[i][13:16].strip() == 'C10'):

            l_nr = 15-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # C4 -> C4
        elif(replace_list[i][13:16].strip() == 'C4'):
            l_nr = 9-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # C3 -> C3
        elif(replace_list[i][13:16].strip() == 'C3'):
            l_nr = 8-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # C2 - C2
        elif(replace_list[i][13:16].strip() == 'C2'):
            l_nr = 6-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # C8 -> C8
        elif(replace_list[i][13:16].strip() == 'C8'):
            l_nr = 13-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # C11 -> C11
        elif(replace_list[i][13:16].strip() == 'C11'):
            l_nr = 16-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # C7 -> C7
        elif(replace_list[i][13:16].strip() == 'C7'):
            l_nr = 12-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # C1 -> C1
        elif(replace_list[i][13:16].strip() == 'C1'):
            l_nr = 2-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # C5 -> C5
        elif(replace_list[i][13:16].strip() == 'C5'):

            l_nr = 10-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # C6 -> C6
        elif(replace_list[i][13:16].strip() == 'C6'):

            l_nr = 11-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H1 -> H2
        elif(replace_list[i][13:16].strip() == 'H1'):
            l_nr = 20-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H2 -> H3
        elif(replace_list[i][13:16].strip() == 'H2'):

            l_nr = 21-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H3 -> H4
        elif(replace_list[i][13:16].strip() == 'H3'):

            l_nr = 22-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H4
        elif(replace_list[i][13:16].strip() == 'H4'):

            l_nr = 23-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line



        # H5 
        elif(replace_list[i][13:16].strip() == 'H5'):

            l_nr = 24-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H6
        elif(replace_list[i][13:16].strip() == 'H6'):

            l_nr = 25-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # H7
        elif(replace_list[i][13:16].strip() == 'H7'):

            l_nr = 26-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H8
        elif(replace_list[i][13:16].strip() == 'H8'):

            l_nr = 27-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # H9
        elif(replace_list[i][13:16].strip() == 'H9'):


            l_nr = 28-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line



        # H10
        elif(replace_list[i][13:16].strip() == 'H10'):

            l_nr = 29-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line
        ########
        # H11
        elif(replace_list[i][13:16].strip() == 'H11'):
            l_nr = 30-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # H12
        elif(replace_list[i][13:16].strip() == 'H12'):
            l_nr = 31-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # H13
        elif(replace_list[i][13:16].strip() == 'H13'):

            l_nr = 32-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # H14
        elif(replace_list[i][13:16].strip() == 'H14'):

            l_nr = 33-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # H15
        elif(replace_list[i][13:16].strip() == 'H15'):

            l_nr = 34-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H16
        elif(replace_list[i][13:16].strip() == 'H16'):

            l_nr = 35-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H17
        elif(replace_list[i][13:16].strip() == 'H17'):

            l_nr = 36-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H18
        elif(replace_list[i][13:16].strip() == 'H18'):

            l_nr = 37-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # H19
        elif(replace_list[i][13:16].strip() == 'H19'):

            l_nr = 38-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H20
        elif(replace_list[i][13:16].strip() == 'H20'):

            l_nr = 39-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H21
        elif(replace_list[i][13:16].strip() == 'H21'):

            l_nr = 40-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        # H22
        elif(replace_list[i][13:16].strip() == 'H22'):

            l_nr = 41-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line


        elif(replace_list[i][13:16].strip() == 'H23'):

            l_nr = 42-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H24
        elif(replace_list[i][13:16].strip() == 'H24'):

            l_nr = 43-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

            
        # H25
        elif(replace_list[i][13:16].strip() == 'H25'):

            l_nr = 44-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H26
        elif(replace_list[i][13:16].strip() == 'H26'):

            l_nr = 45-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line

        # H27
        elif(replace_list[i][13:16].strip() == 'H27'):
            l_nr = 46-1
            line = correct_list[l_nr][0:26]+replace_list[i][26:55]+correct_list[l_nr][55:]

            replace_atomname_residuename[l_nr] = line



    replace_atomname_residuename[16] = correct_list[16]
    replace_atomname_residuename[17] = correct_list[17]
    replace_atomname_residuename[18] = correct_list[18]
    replace_atomname_residuename[46] = correct_list[46]

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
