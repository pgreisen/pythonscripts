from getligandposition import *
from matcherposformat import *
from getconstraintatoms import *


import argparse

"""
Generate position file for design for all residues XX Aangstrom
away from ligand - default is 6

@requires:

python get_positions_around_ligand.py PDBFILE POSITIONFILE LIGANDATOMNAME or CST-file

PDBFILE  : contains coordinates along with the ligand present
POSITIONFILE   : This will be the name of the new position file for each cst
LIGANDATOMNAME  : Ligand atomname e.g., N1
CST : a cst-file with the different constraints 

here is an example of the real cmd

getligandpositionmatchermain.py -f tmp/UM_15_W17Y295_1gca_nohet_1_match_1__DE_1.pdb -p newpos.pos -n O1 N1


@return:
the positions with the cutoff

"""


# @requires: pdbfile, position file, and ligand name
def main():

    parser = argparse.ArgumentParser(description='Narrows down position')

    parser.add_argument("-f", "--file", dest="pdbfile",
                      help="File with aligned ligand")
    
    parser.add_argument("-p", "--pos",
                      dest="pos",help="Position file")

    parser.add_argument("-n", "--name",nargs="*",action="append",
                      dest="ligandnames",help="Ligand names")

    parser.add_argument("-c", "--cons",
                        dest="cst_file",help="Constraint file")

    input_var = parser.parse_args()
    
    gp = GetLigandPositions()
    mp = matcherposformat()
    # Getting the positions 
    #############################positions = gp.get_positions(input_var.pos)
    # Get pdbfile
    p_pdbfile = gp.get_pdbfile(input_var.pdbfile)
    # Get ligand coordinates
    new_positions = []
    if(input_var.ligandnames != None):
        for lig in input_var.ligandnames[0]:
            tmplig = lig.split('-')
            if tmplig[0] == 'aro' and len(tmplig) > 2:
                tmplig.pop(0)
                lig_coor = gp.get_ligandatom_aromatic_coordinates(p_pdbfile,tmplig)
                new_positions.append(gp.get_protein_positions(p_pdbfile,lig_coor))
            else:
                lig_coor = gp.get_ligandatom_coordinates(p_pdbfile,lig)
                new_positions.append(gp.get_protein_positions(p_pdbfile,lig_coor))

    # atom names, number of constratins are taken from the cst-file
    else:
      gc = getconstraintatoms()
      nr,names =  gc.getnumberconstraints(input_var.cst_file)
      for lig in names:
          lig_coor = gp.get_ligandatom_coordinates(p_pdbfile,lig)
          new_positions.append(gp.get_protein_positions(p_pdbfile,lig_coor,positions))


    mp.set_ncst_string(new_positions)


if __name__ == "__main__":
    main()
