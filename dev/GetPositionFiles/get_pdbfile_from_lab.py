#!/usr/bin/env python
import sys, shutil, os, subprocess
import os,sys, argparse
from numpy import *
from numpy import linalg as LA
from collections import defaultdict



class GetPositionFiles:


    def __init__(self):
        self.pdbname = ""
        self.path = ""
        self.cutoff_positionfile = 2


    def set_path_to_files(self):
        self.path = '/lab/shared/scaffolds/'+self.pdbname[1:3]+'/'+self.pdbname+'/'

    def main(self):

        parser = argparse.ArgumentParser(description="Generate Position files for either matching or patchdock")
        # get the initial rosetta design as input
        parser.add_argument("-f", dest="pdbinput", help="PDB ID as input")

        input_variables = parser.parse_args()

        self.pdbname = input_variables.pdbinput

        self.set_path_to_files()

        try:
            if( os.path.exists( self.path ) ):
                # list files
                files = os.listdir( self.path )
                for fl in files:
                    tmpfile = fl.split('.')
                    import pdb; pdb.set_trace()
                    if( tmpfile[0][-1] < self.cutoff_positionfile ):

                        print fl

        except:
            print "The "

        # subprocess.Popen(move_files,shell=True).wait()

if __name__ == "__main__":
    run = GetPositionFiles()
    run.main()
