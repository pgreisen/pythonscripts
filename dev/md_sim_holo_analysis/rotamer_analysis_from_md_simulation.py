from getligandposition import *
from matcherposformat import *
from getconstraintatoms import *
from amber_analysis_format import *
from alignment import *
from amber_ptraj_analysis import *
import argparse, subprocess,os


class RotamerAnalysisfromMDSimulation:


    def __init__(self):
        gp = GetLigandPositions()
        mp = matcherposformat()
        aa = AmberAnalysisFormat()
        self.chi_script = "/Users/pgreisen/pythonscripts/dev/chi_value_analysis"
        self.path_to_amber = "/Users/pgreisen/Programs/amber14/bin"
        self.apo = ""
        self.pdbfile = ""
        self.ligandname = ""

    """
    Generate position file for design for all residues XX Aangstrom
    away from ligand and having an angle less then - default is 6


    @requires:

    # Make an analysis of the holo state of the protein-ligand complex
    python ~/pythonbin/dev/md_sim_holo_analysis/rotamer_analysis_from_md_simulation.py -f Minimization/0BCG_1sjw_20reverted2.pdb -n BCG



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

    def get_first_shell_from_rosetta_design(self,pdbfile,ligandname):
        # NEED TO IMPROVE THIS
        gp = GetLigandPositions()
        mp = matcherposformat()
        aa = AmberAnalysisFormat()

        # Get pdbfile
        p_pdbfile = gp.get_pdbfile_nowater(pdbfile)
        # Get ligand coordinates
        protein_positions = []
        lig_coor = gp.get_ligandatom_coordinates_by_residuename(p_pdbfile,ligandname)

        protein_positions = gp.get_protein_positions(p_pdbfile,lig_coor)

        string_positions = self.write_ligand_positions_to_file( protein_positions )
        return protein_positions


    # @requires: pdbfile, position file, and ligand name
    def main(self):

        print "How to run the script"
        print "python rotamer_analysis_from_md_simulation.py -f UM_1_S137F135T57_1fen_match4_11_.pdb -n SRO -a false"


        parser = argparse.ArgumentParser(description='Analysis and comparison of MD simulation and Rosetta initial design')

        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="pdbfile",
                      help="Initial Rosetta design")


        parser.add_argument("-n", "--name", dest="ligandname",
                      help="Residue name of ligand in Rosetta design")

        # action="store_true" is used envoke the right state of the boolean
        parser.add_argument("-a", "--apo", dest="apo", help="Run analysis on apo if -a apo is given as input")

        input_var = parser.parse_args()

        args_dict = vars(parser.parse_args())
        for item in args_dict:
            setattr(self, item, args_dict[item])



        first_shell_residues_from_rosetta_design = self.get_first_shell_from_rosetta_design(self.pdbfile,self.ligandname)

        ##import pdb; pdb.set_trace()

        positions_in_protein = ""
        for key in first_shell_residues_from_rosetta_design:
            positions_in_protein = positions_in_protein+" "+key

        # chi-angle analysis
        # Just called script and it should dump table
        exe = "python "+self.chi_script+"/chi_analysis.py -f "+self.pdbfile+" -l \" "+positions_in_protein+" \" "
        subprocess.Popen(exe,shell=True).wait()

        # dumping pdb file from MD simulation
        mda = amber_ptraj_analysis()
        mda.get_parameterfile_and_rst_file()
        mda.get_ptraj_analysis_file( self.ligandname, first_shell_residues_from_rosetta_design, self.apo)




if __name__ == "__main__":
   pos = RotamerAnalysisfromMDSimulation()
   pos.main()