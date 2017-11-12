#!/usr/bin/env python
import sys,os, shutil, subprocess, time


'''
The script runs different protocols:

1. Scores the input pdb ( sasa, shape complementarity, IFE, packstat )

2. Number of mutations and which mutations

3. Rotamer Boltzmann for first shell residues

A directory with pdb files - make sure that the naming is BLA_PDBID such that the 
string split one is pdbid[1] - if different set the variable in line 137

python /Users/greisen/pythonbin/dev/post_analysis/post_analysis.py 


'''

# Requires a pdb file
# Returns sequence 
def get_seq(filename):
    tmp = open(filename,'r')
    resid = ''
    seq = []
    j = 0
    for line in tmp:
        tm = line.split()
        if tm[0] == 'ATOM' and j == 0:
            resid = tm[5]
            seq.append(tm[3]+' '+tm[5])
        j = j+1
        if tm[0] == 'ATOM' and tm[5] != resid:
            seq.append(tm[3]+' '+tm[5])
            resid = tm[5]
    return seq


def get_seq_pdbfile(filename):
    tmp = open(filename,'r')
    resid = '0'
    seq = []
    j = 0
    for line in tmp:
        tm = line.split()
        if line[0:4] == 'ATOM' and resid != line[22:26]:
            resid = line[22:26]
            seq.append(line[17:20]+' '+tm[5])
    return seq



# Requires two arrays (same length right now)
# Wild type = a and b = rosetta output
# Prints the differences between them
# Prints number of mutations
# Returns the a list of the mutant's mutations
def print_def(a,b):
    output = open("mutations.dat",'w')
    string = ''
    mt = []
    ln = len(a)
    i = 0
    mut = 0
    output.write("Design     Native\n")
    for i in range(ln):
        if a[i] != b[i]:
            output.write(str(a[i])+"    "+b[i]+"\n")

            mut = mut +1
            mt.append(b[i])

            tmp = a[i].split()
            string = string + '+' + str(tmp[1])

    output.write("For pymol: \n")
    output.write(string)
    return mt, mut


def get_mutations(pdb1,pdb2):

    a = get_seq_pdbfile(pdb1)
    b = get_seq_pdbfile(pdb2)

    seq_a = len(a)
    seq_b = len(b)
    
    tmp = print_def(a,b)


def setup():

    # List the files in the directory
    pdbfiles = os.listdir('./')
    
    dummy = 1

    for pdbfile in pdbfiles:
    
        if pdbfile[-3:] == 'pdb':

            dst = str(dummy)

            os.mkdir(dst)

            shutil.move(pdbfile,dst)

        dummy += 1


def get_differences():

    path = './'

    # List the files in the directory
    dirs = os.listdir(path)
    
    for dir in dirs:
        
        if( os.path.isdir(dir) ):

            os.chdir(dir)

            pdbfiles = os.listdir(path)
        
            for pdbfile in pdbfiles:

                if pdbfile[-3:] == 'pdb':

                    pymol_file = dir+".pdb"
                    
                    # copy files
                    shutil.copy(pdbfile,pymol_file)
                    
                    # Get pdb id
                    pdbid = pdbfile.split('_')[4]

                    pdbid = pdbid.lower()
                    
                    execute = 'python ~/pythonbin/get_pdbfile_from_lab.py '+pdbid


                    subprocess.Popen(execute,shell=True).wait()

                    # native pdb
                    native_pdb = pdbid+"_nohet_1.pdb"

                    # Write a file with the mutations
                    get_mutations(pdbfile,native_pdb)
                
                    # run scores
                    run_score = '~/files/ScoreEnzdes_SC_IFE_PackStat/run.sh '+pdbfile+' > log &'

                    subprocess.Popen(run_score,shell=True).wait()
                    time.sleep(10)
                    # run rotamer Bolztman
                
                    run_rb = '~/files/RotamerBoltzmann/run.sh '+pymol_file+' > rb_log &'
                
                    subprocess.Popen(run_rb,shell=True).wait()
                    time.sleep(10)
            os.chdir('../')


    

def main():

    print "YOU NEED TO SET THE PARAMETER FILES FOR THE SCORE AND THE ROTAMER BOLZTMANN"


    # Setting up a pdb-file into each directory
    setup()
    
    # Get mutations as well as the native pdb
    get_differences()



if __name__ == "__main__":
    main()
