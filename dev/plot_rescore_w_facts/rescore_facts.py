#!/usr/bin/env python
import sys
from pylab import *
import pylab as plt
import numpy as np
import argparse
from collections import defaultdict
import pdb;

class RescoreFacts:


    def __init__(self):
        self.aa_binding_energy = defaultdict(list)
        self.file1 = ""
        self.file2 = ""
        self.scoreterm1 = "ligand_rms_no_super_X"
        self.scoreterm2 = "interface_delta_X"
        self.scoreterm3 = "ddg"
        self.sx = 0
        self.sy = 0
        self.sz = 0
        self.st1 = []
        self.st2 = []
        self.st3 = []
        self.st4 = []

    def set_scores(self):
        with open(self.file1) as f:
         i = 0
         for line in f:
             # Skip the first line
             if( line[0:4] == "SEQU" ):
                 pass
             elif( line[0:4] == "SCOR" and i == 1):
                 indexer = 0
                 x_index = 0
                 y_index = 0
                 tag_index = 0

                 tmpline = line.split()
                 for score in tmpline:
                     if ( score == self.scoreterm1 ):
                         x_index = indexer
                     elif( score == self.scoreterm2 ):
                         y_index = indexer
                     elif( "description" == score ):
                         tag_index = indexer
                     indexer += 1

                 assert tag_index != x_index != y_index

             elif( line[0:4] == "SCOR" and i != 1 ):
                 tmpline =  line.split()
                 self.st1.append( tmpline[ x_index ] )
                 self.st2.append( tmpline[ y_index ] )
                 self.st3.append( tmpline[ tag_index ] )
             i += 1

    def set_rescore(self):
        with open(self.file2) as f:
         i = 0
         for line in f:
             # Skip the first line
             if( line[0:4] == "SEQU" ):
                 pass
             elif( line[0:4] == "SCOR" and i == 1):
                 indexer = 0
                 x_index = 0
                 tmpline = line.split()
                 for score in tmpline:
                     if ( score == self.scoreterm3 ):
                         x_index = indexer

                     indexer += 1

             elif( line[0:4] == "SCOR" and i != 1 ):
                 tmpline =  line.split()
                 self.st4.append( tmpline[ x_index ] )

             i += 1

    def plot_data(self,filename="results"):
        plot(self.st1, self.st2,'.',c='r',markersize=12)
        plot(self.st1, self.st4,'.',c='k',markersize=12)

        # suptitle(name_for_plot)
        # Added for illustration
        #plot(4.085,-15.733,'o',c='r',markersize=14)
        #plot(9.805,-15.134,'o',c='r',markersize=14)
        savefig(filename+".png")
        show()


    def main(self):
        parser = argparse.ArgumentParser(description="Plot results from rescore with facts")
        parser.add_argument('--file1',dest='file1', help='initial file',)
        parser.add_argument('--file2',dest='file2', help='rescored file',)

        #parser.add_argument('--scoreterm1',dest='scoreterm1', help='score term analyzed')
        #parser.add_argument('--scoreterm2',dest='scoreterm2', help='score term analyzed')
        #parser.add_argument('--scoreterm3',dest='scoreterm3', help='score term analyzed')

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.set_scores()
        self.set_rescore()

        self.plot_data()


if __name__ == "__main__":
    run = RescoreFacts()
    run.main()