#!/usr/bin/env python
import sys,shutil
import os

def get_file(filename):
    f = open(filename,'r')
    fl = f.readlines()
    f.close()
    return fl

def read_file(pdbfile,cfactor):
    '''

    Making sure the pdb chain is number from 1 - N 
    without any missing residues. This is important for 
    patchdock and matching.

    '''
    sequence_number = []
    outname = 'tmp.pdb'
    tmp = open(outname,'w')
    # Initialize from 0
    previuos_number = 0

    for line in pdbfile:
        tm = str(line)
        tot = tm
        if tm[0:4] == 'ATOM' and tm[13:15] == 'N ':

            previuos_number = previuos_number + 1
            cf = str( previuos_number )
            if len(cf) == 1:
                cf = '   '+cf
            if len(cf) == 2:
                cf = '  '+cf
            if len(cf) == 3:
                cf = ' '+cf

            first = tm[0:22]
            second = tm[27:]
            tot = first+cf+' '+second
            tmp.write(tot)
        elif tm[0:4] == "ATOM" :
            cf = str( previuos_number )
            if len(cf) == 1:
                cf = '   '+cf
            if len(cf) == 2:
                cf = '  '+cf
            if len(cf) == 3:
                cf = ' '+cf

            first = tm[0:22]
            second = tm[27:]
            tot = first+cf+' '+second
            tmp.write(tot)

        elif( tm[0:4] == 'HETA' ):
            tmp.write(line)

    tmp.close()
    return outname

def get_number_to_substract(pdbfile):

    number = 0

    for line in pdbfile:

        if line[0:4] == 'ATOM':

            number = int(line[23:26]) - 1

            break

    return number


def main():

    filename = sys.argv[1]

    fl = get_file(filename)

    number = get_number_to_substract(fl)

    print "number",number

    if(number == 0):
        print "PDB starts from residue 1"
        src  = read_file(fl,number)
        #pathname = os.path.abspath('')
        #shutil.move(pathname+'/'+str(src),pathname+'/'+str(filename))

    else:
        print "The pdb chain starts at ", number
        src  = read_file(fl,number)
        pathname = os.path.abspath('')
        shutil.move(pathname+'/'+str(src),pathname+'/'+str(filename))
     
if __name__ == "__main__":
    main()
