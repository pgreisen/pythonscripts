import os, shutil,sys
import subprocess
from subprocess import Popen, PIPE, STDOUT

"""

The script sets up directories for the different charge sets and runs am1-bcc computation using AmberTools. 

python setup_electrostatics_benchmark_protocol.py


"""


class Convert_AM1BCC_to_Rosetta():


    def __init__(self):
        self.debug = False


    def convert_amber_atomtype_to_rosetta_atomtype(self):
        """
        @requieres that ligand_am1_bcc.mol2 files is already generated

        @returns a mol2 with atomtypes readable for Rosetta's mol_to_parameter scripts

        """

        tmpfile = open("tmp.mol2", 'w')
        with open("ligand_am1_bcc.mol2",'r') as f:
            atoms = False

            for line in f:
                if( self.debug ):
                    print "ATOM", line.find("@<TRIPOS>ATOM"),line
                    print "BOND", line.find("@<TRIPOS>BOND"),line

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
        self.convert_amber_atomtype_to_rosetta_atomtype()


if __name__ == '__main__':
    run = Convert_AM1BCC_to_Rosetta()
    run.main()