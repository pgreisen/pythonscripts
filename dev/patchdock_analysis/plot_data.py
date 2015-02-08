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
        """
        @requires a list of values
        @returns mean and standard deviations
        """
        return mean(values),sqrt(var(values))

        
    def plot_histogram(self,values):
        hist(values)
        savefig('histogram.png')
        show()


    def plot_histogram_w_mean(self,values,mn,sd):
        hist(values)
        xmin,xmax = xlim()
        ymin,ymax = ylim()
        axvline(mn,c='r', linewidth=4, linestyle='--')
        axvline(mn+sd,c='k', linewidth=4, linestyle='--')
        text( xmin + (xmax-xmin)*0.02,  (ymax-ymin)/2,"Mean = "+str(round(mn,1))+"\n SD= "+str(round(sd,1)),size='x-large')
        savefig('histogram.png')
        show()


    def write_data_to_file(self,values,maxnumber=10):
        maxnumber = len(values)
        filename = open('data.txt','w')
        for val in range(maxnumber):
            filename.write(str(values[val])+'\n')

    def get_top_dictionary(self,dictionary,number=10):
        return dict(dictionary.items()[:number])

    def get_cutoff_value(self,values,cutoff=2):
        # Debug 
        print mean(values)
        print cutoff*sqrt(var(values))
        #############################
        return  mean(values)+cutoff*sqrt(var(values))

    def write_values_above_cutoff(self,dictionary,cutoff,maxpdbs=20):
        datafile = open('get_patchdock_results.sh','w')
        dummy = 0
        for k, v in dictionary.items():
            
            if float(v) >= cutoff and dummy <= maxpdbs:
                pdbname,N = k.split('_')

                # perl /work/greisen/ExternalProgram/patchdock/PatchDock/transOutput.pl $1 $2 $3;

                template = 'perl /work/greisen/ExternalProgram/patchdock/PatchDock/transOutput.pl '+str(pdbname)+'_nohet_1_patchdock.out '+str(N)+' '+str(N)+' #  '+str(v)+'\n'

                datafile.write(template)
                dummy = dummy + 1

    def write_table_to_file(self, hashtable):
        """ write the content of hash table to file"""
        filename = open('table.dat','w')
        for key in hashtable:
            filename.write(str(key)+","+str(hashtable[key])+"\n")

