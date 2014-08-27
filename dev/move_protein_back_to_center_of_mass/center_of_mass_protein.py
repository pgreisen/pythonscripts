#!/usr/bin/python
from optparse import OptionParser
from math import sqrt
from numpy import *
import os,shutil

'''

'''

# Requires PDB files
# Returns list with lines in file
def get_pdbfile(pdbfile):
    fl = open(pdbfile,'r')
    pdb = fl.readlines()
    fl.close()
    return pdb

def get_center_of_mass(chainid, pdbfile ):
    center_of_mass = array([0.0,0.0,0.0])
    number_of_atoms = 0

    for line in pdbfile:
        if ( line[0:4] == 'ATOM' and str(line[20:22]).strip() == chainid ):
            center_of_mass += array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
            number_of_atoms += 1
    com = center_of_mass/number_of_atoms
    print "The center of mass is ", com
    return com

# Requires a number
# Set the right length for the number
def set_length_digit(number):
    number = str(number)
    lngth = len(number)
    if lngth == 7:
        return ' '+number
    if lngth == 6:
        return '  '+number
    if lngth == 5:
        return '   '+number
    if lngth == 4:
        return '    '+number
    else:
        return number


def move_pdb_to_center_of_mass(com,pdbfile,chainid):
    with open("com_"+chainid+".pdb",'w') as f:
        for line in pdbfile:
            if ( line[0:4] == 'ATOM' and str(line[20:22]).strip() == chainid ):
                tmp = array([float(line[31:39]),float(line[39:47]),float(line[47:55])]) - com
                x = '%.2f' %tmp[0]
                y = '%.2f' %tmp[1]
                z = '%.2f' %tmp[2]
                x =  set_length_digit(x)
                y =  set_length_digit(y)
                z =  set_length_digit(z)
                newline = line[0:30]+str(x)+str(y)+str(z)+line[55:]
                f.write(newline)

def main():

    parser = OptionParser()

    parser.add_option('-f',dest='PDB',
                      help='Cleaned pdb file with metal ion present default=ZN')
    parser.add_option('-c',dest='chainid',
                      help='Chain id for protein chain default=A', default='A')



    (options,args) = parser.parse_args()
    
    pdbfile = get_pdbfile(options.PDB)

    chainid = options.chainid

    com = get_center_of_mass(chainid, pdbfile )

    move_pdb_to_center_of_mass(com,pdbfile,chainid)



if __name__ == "__main__":
    main()
