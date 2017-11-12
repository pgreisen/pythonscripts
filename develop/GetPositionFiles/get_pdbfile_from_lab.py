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

        self.positions = []



    def set_path_to_files(self):
        self.path = '/lab/shared/scaffolds/'+self.pdbname[1:3]+'/'+self.pdbname+'/'


    def set_positions(self,positionfile):
        # print positionfile
        with open(self.path+positionfile,'r') as f:
            for line in f:
                tmp = line.split()
                for pos in tmp:
                    self.positions.append(pos)


    def write_positionfile(self, postions_in_set):
        with open("pos.pos",'w') as f:
            for pos in postions_in_set:
                f.write(pos+" ")




    def main(self):

        parser = argparse.ArgumentParser(description="Generate Position files for either matching or patchdock")
        # get the initial rosetta design as input
        parser.add_argument("-f", dest="pdbinput", help="PDB ID as input")

        input_variables = parser.parse_args()

        self.pdbname = input_variables.pdbinput

        self.set_path_to_files()

        #try:
        if( os.path.exists( self.path ) ):
            # list files
            files = os.listdir( self.path )
            for fl in files:
                # print fl
                tmpfile = fl.split('.')
                #import pdb; pdb.set_trace()
                if(  fl.endswith(".pos") and int( tmpfile[0][-1] ) < self.cutoff_positionfile ):
                    # print fl
                    self.set_positions( fl )



        #except:
        #    print "The file didnt exists"
        self.write_positionfile( set(self.positions) )
        # subprocess.Popen(move_files,shell=True).wait()

if __name__ == "__main__":
    run = GetPositionFiles()
    run.main()
