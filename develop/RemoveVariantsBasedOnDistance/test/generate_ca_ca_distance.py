#!/usr/bin/env python
import sys
from numpy import *
import argparse

"""

"""


class CACA_distance:


    def __init__(self):
        self.cutoff = 8.0
        self.pdbfile = ""
        self.ca_coordinates = {}
        self.caca_distance = {}
        self.cutoffatom = "CB"
        self.ca_coor = {}
        self.aa = {
            'A' : "ALA",
            'R' : "ARG",
            'N' : "ASN",
            'D' : "ASP",
            'C' : "CYS",
            'Q' : "GLN",
            'E' : "GLU",
            'G' : "GLY",
            'H' : "HIS",
            'I' : "ILE",
            'L' : "LEU",
            'K' : "LYS",
            'M' : "MET",
            'F' : "PHE",
            'P' : "PRO",
            'S' : "SER",
            'T' : "THR",
            'W' : "TRP",
            'Y' : "TYR",
            'V' : "VAL",
            '-' : "DEL"
        }

        self.aa3letter = {
             "ALA" : 'A',
             "ARG" : 'R',
             "ASN" : 'N',
             "ASP" : 'D',
             "CYS" : 'C',
             "GLN" : 'Q',
             "GLU" : 'E',
             "GLY" : 'G',
             "HIS" : 'H',
             "ILE" : 'I',
             "LEU" : 'L',
             "LYS" : 'K',
             "MET" : 'M',
             "PHE" : 'F',
             "PRO" : 'P',
             "SER" : 'S',
             "THR" : 'T',
             "TRP" : 'W',
             "TYR" : 'Y',
             "VAL" : 'V',
             "DEL" : '-' 
        }

        self.remove = []
        self.debug = 0

    # return coordinates of residue
    def set_calpha_coordinates( self ):
        with open(self.pdbfile,'r') as f:
            for line in f:
                
                if( line[0:4] == "ATOM" and line[17:20] == "GLY" and line[13:15] == "CA"):
                    print "We have a glycine residue in the list"
                    key = 'G'+line[23:26].strip()+'_'+line[21:22].strip()
                    x = str(line[30:38]).rstrip()
                    y = str(line[38:46]).rstrip()
                    z = str(line[46:54]).rstrip()
                    self.ca_coor[key] =  array([float(x),float(y),float(z)])
                elif( line[0:4] == "ATOM" and line[13:15] == self.cutoffatom):
                    key = self.aa3letter[ line[17:20] ]+str(line[23:26].strip())+'_'+line[21:22].strip()
                    x = str(line[30:38]).rstrip()
                    y = str(line[38:46]).rstrip()
                    z = str(line[46:54]).rstrip()
                    self.ca_coor[key] =  array([float(x),float(y),float(z)])


    def set_caca_distance(self):

        for key in self.ca_coor:
            for key2 in self.ca_coor:

                tmp_length = linalg.norm( self.ca_coor[key] - self.ca_coor[key2])
                tmpkey = key+","+key2

                if( tmp_length < self.cutoff and key != key2):
                    tmptmp = tmpkey+","+str(tmp_length)
                    self.remove.append(tmptmp)
                self.caca_distance[tmpkey] = tmp_length

    def set_keys(self, datafile):
        with open(datafile, 'r') as f:
            for line in f:
                tmpline = line.split('_')
                tmp_res = tmpline[0][0]
                tmp_position = tmpline[0][1:-1]
                tmp_chain = tmpline[1][0]
                tmp_key = tmp_res+tmp_position+"_"+tmp_chain
                tmp_value = (self.aa[tmp_res], tmp_position, tmp_chain )
                self.ca_coordinates[tmp_key] = tmp_value

        if(self.debug == 1):
            print "The following ca-ca coordinates were read ", self.ca_coordinates

    def write_to_file(self):
        with open("distance_matrix",'w') as f:
            for key in self.caca_distance.keys():
                f.write(key+","+str( round( self.caca_distance[key] ,3))+"\n")

    def main(self):


        # File name of pdb file
        #self.datafile = sys.argv[1]
        self.pdbfile = sys.argv[1]
        #self.set_keys(self.datafile)
        self.set_calpha_coordinates()
        self.set_caca_distance()
        self.write_to_file()

        if(self.debug == 1):
            print "The collected ca-coordinates: ", self.ca_coordinates


        with open("cbcb_distance_cutoff.dat",'w') as f:
            for i in self.remove:
                f.write(i+"\n")


        '''
        x = str(line[30:38]).rstrip()
        y = str(line[38:46]).rstrip()
        z = str(line[46:54]).rstrip()
        tmp_vector = array([float(x),float(y),float(z)])
        tmp_length = linalg.norm(tmp_vector - i)
        '''


if __name__ == "__main__":
   run = CACA_distance()
   run.main()
