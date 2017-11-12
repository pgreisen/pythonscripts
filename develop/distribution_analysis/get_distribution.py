
#!/usr/bin/env python
import os,shutil
from optparse import OptionParser
from pylab import *
from numpy import mean,median

def read_file(filename):
    f = open(filename,'r')

    number = [float(number.strip('\n')) for number in f.readlines() ]

    f.close()

    return number

def concatenate_lists(a,b):
    c = a+b
    return c


def main():

    parser = OptionParser()

    parser.add_option('-f',dest='datafile',
                      help='Datafile')

    parser.add_option('-l',dest='list_of_datafiles',
                      help='List containing naming of data files')

    parser.add_option('-c',dest='concatenate',
                      help='Should all the data be pooled together')



    (options,args) = parser.parse_args()


    # String of atom names from the ligand to be aligned with
    # Set the path to the pdb file
    datafile = options.datafile

    list_of_datafiles = options.list_of_datafiles
    
    total_data = []


    concatenate = False

    if(datafile != None):

        datafile = read_file(datafile)

        hist(datafile,100,normed=1)
        
        savefig('histogram.png')


    # Multiple disitrbution plot 

    elif(list_of_datafiles != None):

        fl_list = open(list_of_datafiles,'r')

        f = [f.strip('\n') for f in fl_list.readlines() ]

        fl_list.close()


        if(concatenate == True):
            for ln in f:
                
                datafile = read_file(ln)
                
                total_data = concatenate_lists(total_data,datafile)

            hist(total_data,100,normed=1)

            savefig('histogram.png')
    

        else:
            dummy = 0
            dataset = []

            for ln in f:

                dataset.append(read_file(ln))
                
                print "The mean is ", round(mean(dataset[dummy]),3)
                print "The median is ",round(median(dataset[dummy]),3)

                hist(dataset[dummy],100,normed=1)

                dummy = dummy + 1
    
            savefig('mulitple_histogram.png')
    

if __name__ == "__main__":
    main()
