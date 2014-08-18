#!/usr/bin/env python
import os,shutil
from optparse import OptionParser
from pylab import *
from analyse_patchdock_score import *
from plot_data import *

"""
Makes an alignment to a zinc site in a protein. Atom names of
ligand has to be set in main in the list SUBSTRATEATOMS

@param -f: pdbfile with zincprotein
@type f  : string
@param -a: atom name to be analysed e.g., NE1 in TRP
@type a  : string
@param -t: atomtype which will be the target of the analysis e.g., carbon-atom 
@type a  : string
@param -r: residue which will be analysed 3-letter e.g., TRP 
@type r  : string
@param -m: minimum distance 
@type r  : float
@param -c: maximum distance
@type r  : float

@return

"""

def main():
    parser = OptionParser()

    parser.add_option('-f',dest='pdbfile',
                      help='Protein file')
    parser.add_option('-l',dest='file_pdbfile',
                      help='File with pdbs')

    (options,args) = parser.parse_args()


    # String of atom names from the ligand to be aligned with
    # Set the path to the pdb file
    pdbfile = options.pdbfile

    file_pdbs = options.file_pdbfile
    
    # Initialize class
    ap = analyse_patchdock()

    if(pdbfile != None):

        pdbname = pdbfile.split('_')[2]

        datafile = ap.get_file(pdbfile)

        ap.get_scores(datafile,pdbname)


    elif(file_pdbs != None):

        fl_list = open(file_pdbs,'r')

        f = [f.strip('\n') for f in fl_list.readlines() ]

        fl_list.close()

        for ln in f:

            pdbname =  ln.split('_')[2]  

            pdbfile = ap.get_file(ln)  

            ap.get_scores(pdbfile,pdbname)
            ap.print_number_scores()

    test = ap.get_sorted_results()

    plt = plot_data()
    values = plt.get_values(test)
    plt.write_data_to_file(values)

    cutoff_value =  plt.get_cutoff_value(values)
    print cutoff_value
    cutoff_value = 3431
    plt.write_values_above_cutoff(test,cutoff_value)

    plt.plot_histogram(values)
    
    
if __name__ == "__main__":
    main()
