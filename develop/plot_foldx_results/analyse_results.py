import sys
from pylab import *


def get_single_amino_acids(residuename):
    aa = {
        "ALA" : "A", 
        "ILE" : "I",
        "LEU" : "L",
        "VAL" : "V",
        "MET" : "M",
        "PHE" : "F",
        "TYR" : "Y",
        "ARG" : "R",
        "LYS" : "K",
        "TRP" : "W",
        "ASN" : "N",
        "GLN" : "Q",
        "ASP" : "D",
        "GLU" : "E",
        "SER" : "S",
        "THR" : "T",
        "GLY" : "G",
        "PRO" : "P",
        "CYS" : "C",
        "HIS" : "H"
        }
    return aa[residuename]


def get_datafile(inputfile):
    data = {}
    energies = []
    with open(inputfile,'r') as f:
        for line in f:
            key,value = line.split()
            # print key
            if ( float(value) > 1.0 ):
                data[key] = value
            energies.append(float( value))
    return data, energies



# print get_single_amino_acids("ALA")
datafile = sys.argv[1]
data, energies = get_datafile( datafile )

#hist(energies)
#show()

for i in data:
    print i,data[i]
