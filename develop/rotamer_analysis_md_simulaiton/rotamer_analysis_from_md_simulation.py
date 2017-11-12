from getligandposition import *
from matcherposformat import *
from getconstraintatoms import *
from amber_analysis_format import *
from alignment import *


import argparse, subprocess,os


class RotamerAnalysisfromMDSimulation:


    def __init__(self):
        pass


    """
    Generate position file for design for all residues XX Aangstrom
    away from ligand and having an angle less then - default is 6


    @requires:

    python get_positions_around_ligand.py PDBFILE POSITIONFILE LIGANDATOMNAME or CST-file

    PDBFILE  : contains coordinates along with the ligand present
    POSITIONFILE   : This will be the name of the new position file for each cst
    LIGANDATOMNAME  : Ligand atomname e.g., N1
    CST : a cst-file with the different constraints

    here is an example of the real cmd

    getligandpositionmatchermain.py -f tmp/UM_15_W17Y295_1gca_nohet_1_match_1__DE_1.pdb -p newpos.pos -n O1 N1

    or executed with center of aromatic ring

    @return: file with position for each cst/rst with the name specified by -p

    Positions with the cutoff and an angle pointing towards the ligand

    """

    def write_ligand_positions_to_file(self,protein_positions):
        tmp_string = ''
        tmp_file = open("first_shell_residues.txt",'w')
        #print "The protein positions are ",protein_positions
        #print "The protein positions are ",protein_positions
        for i in protein_positions:
            tmp_string += i+" "
            tmp_file.write(i+" ")
        tmp_file.close()
        return tmp_string

    def write_rmsd_analysis(self,list_with_text):
        outputfile = open("rmsd_of_ligating_residues.dat",'w')
        for i in list_with_text:
            outputfile.write(i+"\n")


    # @requires: pdbfile, position file, and ligand name
    def main(self):

        parser = argparse.ArgumentParser(description='Narrows down position')

        parser.add_argument("-f", "--file", dest="pdbfile",
                      help="File with aligned ligand")

        parser.add_argument("-f2", "--file2", dest="pdbfile2",
                      help="File with aligned ligand")


        input_var = parser.parse_args()
    
        gp = GetLigandPositions()
        mp = matcherposformat()

        aa = AmberAnalysisFormat()


        # Get pdbfile
        # p_pdbfile = gp.get_pdbfile(input_var.pdbfile)


        # Get pdbfile without water present
        p_pdbfile = gp.get_pdbfile_nowater(input_var.pdbfile)

        # Get ligand coordinates
        protein_positions = []

        lig_coor = gp.get_ligandatom_coordinates_by_residuename(p_pdbfile,"750")

        # protein_positions.append(gp.get_protein_positions(p_pdbfile,lig_coor))

        protein_positions = gp.get_protein_positions(p_pdbfile,lig_coor)


        # what to do with the protein positions
        # print "Here we go", protein_positions

        string_positions = self.write_ligand_positions_to_file( protein_positions )

        aa.write_amber_analysis_to_file(protein_positions,p_pdbfile)


        # Need to align the two protein for the analysis
        # One could have small drifts during MD or relax with Rosetta
        al = Alignment()
        al.align_proteinA_onto_proteinB(p_pdbfile, p_pdbfile)


        pdbfile2 = gp.get_pdbfile(input_var.pdbfile2)


        rmsd_analysis = []

        for i in protein_positions:

            tmp = gp.get_coordinates_residue_number_not_hydrogen( p_pdbfile, i )
            tmp2 = gp.get_coordinates_residue_number_not_hydrogen( pdbfile2, i )

            rmsd = al.rmsd(tmp, tmp2)

            tmp_string = protein_positions[i]+" "+str(i)+" "+str(round(rmsd,3))

            rmsd_analysis.append( tmp_string )

        self.write_rmsd_analysis( rmsd_analysis )


        # chi-angle analysis
        # Just called script and it should dump table
        exe = "python ~/pythonbin/chi_analysis.py -f "+input_var.pdbfile+" -l \" "+string_positions+" \" "
        subprocess.Popen(exe,shell=True).wait()

         # chi-angle analysis
        # Just called script and it should dump table
        exe = "python ~/pythonbin/chi_analysis.py -n chi_table_"+input_var.pdbfile2+".dat -f "+input_var.pdbfile2+" -l \" "+string_positions+" \" "
        subprocess.Popen(exe,shell=True).wait()




if __name__ == "__main__":
   pos = RotamerAnalysisfromMDSimulation()
   pos.main()