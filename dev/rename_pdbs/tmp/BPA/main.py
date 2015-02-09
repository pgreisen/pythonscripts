__author__ = 'greisen'

import os,sys

'''

renaming beta-propeller proteins

'''


import argparse, subprocess,os
import os, shutil



# @requires: pdbfile, position file, and ligand name
def main():

    pdbid = ""

    # get the file in the directory
    path = './'
    files = os.listdir(path)

    pdbid = sys.argv[1]

    for pdb in files:
        if(pdb.startswith("patchdock")):
            new_pdbid = pdbid+"_"+pdb
            shutil.move(pdb, new_pdbid)



if __name__ == "__main__":
   main()

