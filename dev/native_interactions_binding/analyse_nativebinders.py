import os,shutil,sys
import pylab as pl


'''
Plots the patchdock score versus the interface energy after rosetta design. 

The script requires that you are in the main directory and each subdirectory contains a 
pdb-file and rosetta output. 

It generates a scatter plot of the data and it is executed like this

python2.7 analyse_patchdock_score.py

'''

dic_sc = {}
dic_ife = {}
dic_packstat = {}
dic_dsasa = {}
dic_ddg = {}



# Requires scorefile 
# Returns list with lines
def get_file(file):
    fl = open(file,'r')
    pfile = fl.readlines()
    fl.close()
    return pfile



def get_results(file,dir):

    for line in file:

        if( len(line) > 5 and line[0:2] == 'SC'):
            tmp = line.split()
            dic_sc[dir] = tmp[1]

        if( len(line) > 5 and line[0:3] == 'ddg'):
            tmp = line.split()
            dic_ddg[dir] = tmp[1]
            

        if( len(line) > 5 and line[0:3] == 'dsa'):
            tmp = line.split()
            dic_dsasa[dir] = tmp[1]
            
        if( len(line) > 5 and line[0:3] == 'int'):
            tmp = line.split()
            dic_ife[dir] = tmp[1]

        if( len(line) > 5 and line[0:3] == 'pst'):
            tmp = line.split()
            dic_packstat[dir] = tmp[1]



def main():


    data_file = "pdb_chainA_0001.pdb"

    path = './'
    
    dirs = os.listdir(path)
    
    for dir in dirs:

        if os.path.isdir(dir):

            os.chdir(dir)

            subdirs = os.listdir(path)


            for subdir in subdirs:

                if(subdir == "BenchmarkEnzdes"):

                    os.chdir(subdir)

                    files = [f for f in os.listdir('.') if os.path.isfile(f)]

                    for file in files:

                        if file == "pdb_chainA_0001.pdb":

                            tmp_file = get_file(data_file)

                            get_results(tmp_file,dir)

                    os.chdir('../')
            os.chdir('../')
                    
    print dic_sc
    print dic_sc["1A28"]


if __name__ == "__main__":
    main()    
