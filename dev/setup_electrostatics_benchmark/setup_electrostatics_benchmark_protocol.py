import os, shutil,sys
import subprocess
from subprocess import Popen, PIPE, STDOUT

"""

The script sets up directories for the different charge sets and runs am1-bcc computation using AmberTools. 

python setup_electrostatics_benchmark_protocol.py


"""


class SetupElectrostaticsBenchmarkProtocol():


    def __init__(self):
        self.rosetta_charge = "parameters_rosetta"
        self.gasteiger_charge = "parameters_gasteiger"
        self.am1bcc = "parameters_am1bcc"

    def get_missing_files(self,dst):
        print "get_missing_files"

        file_is_present = False

        # self.rosetta_charge
        # cd parameters_rosetta
        os.chdir( self.rosetta_charge )

        # loop over all the files
        path = './'
        files = os.listdir( path )
        for fl in files:
            if(fl == "ligand.mol2"):
                file_is_present = True
                break

        if file_is_present == False :
            exe = "~nivon/install/openbabel/bin/babel -imol ligand.mol -omol2 ligand.mol2"
            p = subprocess.Popen(exe,stdin=PIPE, stdout=PIPE, stderr=STDOUT,shell=True)# .wait()
            file_is_present = True
            stdout, stderr = p.communicate()
            print p, stderr

        os.chdir("../")

        try:
            shutil.copy(self.rosetta_charge+"/ligand.mol2", dst)
        except:
            print "#########################################"
            print "The ligand.mol2 is not generated"
            sys.exit()



    def setup_am1_bcc_and_generate_rosetta_parameters(self):
        try:
            # mkdir directories
            os.mkdir(self.am1bcc)
        except:
            print "The directory already exists ", self.am1bcc

        self.get_missing_files( self.am1bcc )

        os.chdir( self.am1bcc )

        exe = "~/amber12/bin/antechamber -i ligand.mol2 -fi mol2 -o ligand_am1_bcc.mol2 -fo mol2 -c bcc"

        # run the AmberTools implementation of the cmd.
        subprocess.Popen(exe,shell=True).wait()

        # change atomtypes from Amber to Rosetta
        self.convert_amber_atomtype_to_rosetta_atomtype()

        # copy to rename file over to the original name 
        shutil.move("tmp.mol2","ligand_am1_bcc.mol2")

        exe_rosetta = "python ~/git_rosetta/Rosetta/main/source/src/python/apps/public/molfile_to_params.py ligand_am1_bcc.mol2 -n LG1 -c --clobber"

        subprocess.Popen(exe_rosetta,shell=True).wait()

        os.chdir("../")


    def convert_amber_atomtype_to_rosetta_atomtype(self):
        """
        @requieres that ligand_am1_bcc.mol2 files is already generated

        @returns a mol2 with atomtypes readable for Rosetta's mol_to_parameter scripts

        """

        tmpfile = open("tmp.mol2", 'w')
        with open("ligand_am1_bcc.mol2",'r') as f:
            atoms = False
            for line in f:

                if ( len(line) > 13 and line.find("@<TRIPOS>ATOM") >-1.0):
                    atoms = True
                elif ( len(line) > 13 and line.find("@<TRIPOS>BOND") >-1.0):
                    atoms = False

                elif( atoms == True and len(line) > 75 ):
                    tmp_characters = line[47]+"."+line[48]
                    line = line[0:47]+tmp_characters+line[50:]

                tmpfile.write(line)
        tmpfile.close()



    def main(self):
        path = './'
        dirs = os.listdir( path )

        # loop over all the directories
        for dir in dirs:
            if ( os.path.isdir ( dir ) and not dir.startswith('.') ):
                # chdir to PDBID NAME
                os.chdir( dir )
                sdirs = os.listdir( path )

                # make sure that gasteiger and rosetta default parameters are
                # present
                subdirs = {}

                for sdir in sdirs:

                    if(os.path.isdir( sdir )):
                        subdirs[ sdir ] = True

                # Setup and run the am1-bcc in AmberTool
                self.setup_am1_bcc_and_generate_rosetta_parameters()

                # make sure the pdb file has the pdbid present in its name
                ##if( subdirs[ self.gasteiger_charge] and subdirs[ self.rosetta_charge] ):
                ##    continue
                ##else:
                ##    print "Directories are missing"
                ##    # sys.exit()

                os.chdir("../")
                print "The path is now ",os.getcwd()


if __name__ == '__main__':
    run = SetupElectrostaticsBenchmarkProtocol()
    run.main()
