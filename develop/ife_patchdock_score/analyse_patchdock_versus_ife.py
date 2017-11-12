import os,shutil,sys
import pylab as pl


'''
Plots the patchdock score versus the interface energy after rosetta design. 

The script requires that you are in the main directory and each subdirectory contains a 
pdb-file and rosetta output. 

It generates a scatter plot of the data and it is executed like this

python2.7 analyse_patchdock_score.py

'''


# Requires scorefile 
# Returns list with lines
def get_file(file):
    fl = open(file,'r')
    pfile = fl.readlines()
    fl.close()
    return pfile


# Requires list, maximum_number
# Returns list with scores in that interval
def get_scores(pfile,numbers):
    pd_scores = []
    start = False
    j =  0
    for line in pfile:
        tmp_line = line.split()
        
        if(len(tmp_line) > 0 and  str(tmp_line[0]) == 'Best'):
            break

        if(start == True):
            if(tmp_line[0] in numbers):
                j = j + 1

                pd_scores.append(tmp_line[2])

        if(len(tmp_line) > 0 and  str(tmp_line[0]) == '#'):
            start = True

    return pd_scores



def get_interface_energy(file):
    f = get_file(file)
    interface_energy = 0
    for line in f:
        tmp_line = line.split()
        if len(tmp_line) > 1 and tmp_line[0] == 'interfE':
            interface_energy = tmp_line[1]

    return interface_energy






# Returns the PD designs
def get_designs_from_pd():
    path = './'
    files = os.listdir(path)
    numbers = []
    interface_energy = []
    pd_design_number = []


    for fl in files:
        
        if os.path.isfile(fl):

            
            tmp =  fl.split('_')

            if( len(tmp) == 6 and fl[-3:] == 'out'):

                datafile = get_file(fl)


            elif( len(tmp) > 6 and fl[-3:] == 'pdb'):
                # lstrip('0')
                index = fl.index('_in_') + 4

                pd_design_number.append( int( fl[index:index+4].lstrip('0') ) )
                
                interface_energy.append(get_interface_energy(fl))

                run = False

                integer = ''
                
                tmp_number = tmp[6].split('.')[0]

                for i in list(tmp_number):

                    if i != '0':

                        run = True
                        
                    if(run):

                        integer = integer + i

                numbers.append(integer)

    pd_scores = get_scores(datafile,numbers)

    return interface_energy, pd_scores, pd_design_number



def main():

    total_interface_energy = []
    total_pd_scores = []
    total_design_numbers = []

    path = './'
    
    dirs = os.listdir(path)
    
    for dir in dirs:

        if os.path.isdir(dir):

            print dir

            os.chdir(dir)

            interface_energy, pd_scores, design_numbers = get_designs_from_pd()

            total_interface_energy = total_interface_energy + interface_energy 

            total_pd_scores = total_pd_scores + pd_scores

            total_design_numbers = total_design_numbers + design_numbers

            assert(len(total_pd_scores) == len(total_interface_energy))

            os.chdir('../')

    #    print total_pd_scores, total_interface_energy
    pl.scatter(total_pd_scores, total_interface_energy)
    pl.savefig('scatterplot.png')
    pl.clf()
    pl.hist(total_design_numbers,bins=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),normed=1)
    pl.xlim(1,15)
    pl.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))

    pl.savefig('design_number.png')
    pl.clf()
    pl.scatter(total_design_numbers,total_interface_energy)
    pl.savefig('design_number_versus_ife.png')

if __name__ == "__main__":
    main()    
