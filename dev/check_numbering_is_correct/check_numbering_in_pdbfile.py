#!/usr/bin/env python
import sys,shutil
import os,subprocess

def get_file(filename):
    f = open(filename,'r')
    fl = f.readlines()
    f.close()
    return fl

def get_pdb_sequence(pdbfile):
    '''

    Making sure the pdb chain is number from 1 - N 
    without any missing residues. This is important for 
    patchdock and matching.

    '''

    sequence = []


    for line in pdbfile:

        tm = str(line)

        if tm[0:4] == 'ATOM':

            cf= int(tm[23:26])

        sequence.append(cf)

    return set(sequence)


def determine_consecutive_sequence(sequence):

    not_consecutive = True

    dummy = 1

    for i in sequence:
        #
        if i != dummy:
            # print i,dummy
            not_consecutive = False
        dummy += 1
    return not_consecutive
        


def main():

    pdbfiles_missing_density = open('missing_density.txt','w')


    pdbfiles_remove_missing_density = open('remove_missing_density.sh','w')
    pdbfiles_remove_missing_density.write('rm -f ')


    dst = 'Fix_sequences'
    try:
        os.mkdir(dst)
    except:
        print "Directory exists"

    path = './'

    dirs = os.listdir(path)

    correct_seq = True

    for dir in dirs:
        
        if os.path.isdir(dir):

            os.chdir(dir)

            filename = 'protein_0001.pdb'

            fl = get_file(filename)

            seq  = get_pdb_sequence(fl)

            correct_seq = determine_consecutive_sequence(seq)

            os.chdir('../')
    
            if(correct_seq == False):
                
                pdbfiles_missing_density.write(dir+' ')

                pdbfiles_remove_missing_density.write(' *'+dir+'* ')

                mvfl = 'mv '+str(dir)+' '+str(dst)+'/'
                
                subprocess.Popen(mvfl,shell=True).wait()
 



if __name__ == "__main__":
    main()
