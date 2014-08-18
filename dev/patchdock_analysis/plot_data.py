#!/usr/bin/env python
from pylab import *
from numpy import *


class plot_data:
    
    # Requires dictionary 
    # @returns numpu.array with values
    def get_values(self,dictionary):
        values = []
        for k, v in dictionary.items():
            values.append(float(v))
        return values


    def get_mean(self,values):
        return mean(values)

        
    def plot_histogram(self,values):
        hist(values)
        savefig('histogram.png')
        show()

    def write_data_to_file(self,values,maxnumber=10):
        maxnumber = len(values)
        filename = open('data.txt','w')
        for val in range(maxnumber):
            filename.write(str(values[val])+'\n')

    def get_top_dictionary(self,dictionary,number=10):
        return dict(dictionary.items()[:number])

    def get_cutoff_value(self,values):
        # Debug 
        print mean(values)
        print 2*sqrt(var(values))
        #############################
        return  mean(values)+2*sqrt(var(values))

    def write_values_above_cutoff(self,dictionary,cutoff,maxpdbs=20):
        datafile = open('get_patchdock_results.sh','w')
        dummy = 0
        for k, v in dictionary.items():
            
            if float(v) >= cutoff and dummy <= maxpdbs:
                pdbname,N = k.split('_')
                template = 'perl ~/files/RetrievePatchdock/transOutput.pl ART_0001.fa_'+str(pdbname)+'_nohet_1_patchdock.out '+str(N)+' '+str(N)+' #  '+str(v)+'\n'

                datafile.write(template)
                dummy = dummy + 1
