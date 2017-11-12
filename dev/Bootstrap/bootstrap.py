#!/usr/bin/env python

import sys
from numpy import *
import argparse
from random import randint
import pylab as plt
from pylab import *

'''
 Expect csv-file with sequence and value:
 EADFEFDA,1.2

'''

class BootstrapData:

    def __init__(self):
        self.nativeseq = ""
        self.seq_data = {}
        self.datafile = ""
        self.native_datafile = ""
        self.bootstrap = 1000
        self.nativesequencefile = ""
        self.variants_sequence = {}
        self.bootstrap_results = {}
        self.offset = 0
        # self.title = "TGT Assay"
        self.title = ""
        #self.title = "F8 Inhibition"
        # self.ylabel = "Counts"
        self.ylabel = ""
        #self.ylabel = "F8 Interference"
        self.mean_reference = 0
        self.chain_name = "LC"
        self.assay = "assay1"
        self.round_nr = 2


    def set_native_sequence(self):
        '''
        set the native sequence used to identify positions that are substituted
        :return:
        '''
        with open(self.nativesequencefile,'r') as f:
            for line in f:
                if (line[0] == ">"):
                    continue
                else:
                    self.nativeseq += line.strip()

    def set_variants_data(self):
        '''
        set hash with variant sequence and data values
        :return:

        2017-01-09
        list with values as redundancy among heavy/light chains
        e.g. same light chain with different heavy chains etc.

        '''
        with open(self.datafile,'r') as f:
            for line in f:
                if(line[0] == ">"):
                    break
                tmpline = line.strip().split(',')
                # key is the sequence
                key = tmpline[0]
                value =  tmpline[1]

                if(key not in self.seq_data.keys() ):
                    self.seq_data[key] = []
                self.seq_data[key].append( value )


    def set_difference(self):
        '''

        assume not deletions or insertions. Set a hash table with positions
        as first key and amino acid as second key with the sequence as a value

        '''
        seqlength = len(self.nativeseq)
        # adjust counting
        dummy = 1 + self.offset
        mutations = []
        for i in range(seqlength):
            # loop over sequences as keys
            for key in self.seq_data:
                # print len(self.nativeseq), len(key),i
                if (self.nativeseq[i] != key[i]):
                    newkey = str(dummy)+"_"+key[i]

                    # position is added as a key to bootstrap
                    if( newkey not in self.variants_sequence.keys() ):
                        self.variants_sequence[newkey] = []
                    # AA at position is added as key with sequence as value
                    #if( key[i] not in self.variants_sequence[str(dummy)].keys() ):
                    #    self.variants_sequence[str(dummy)][key[i]] = []
                    # is this the smartest way of doing it
                    self.variants_sequence[newkey].append(key)
            dummy += 1


    def bootstrap_values(self):
        '''

        sets bootstrapped values for each of the positions
        and amino acids (self.bootstrap_results)


        '''
        # self.bootstrap_results
        # loop over residue number and aa substitution
        for residuenr in self.variants_sequence.keys():
            # loop over substitution at this position
            nr_variants = len(self.variants_sequence[residuenr])

            if(nr_variants > 1 ):
                newkey = residuenr
                if (newkey not in self.bootstrap_results.keys()):
                    self.bootstrap_results[newkey] = []

                # start bootstrapping
                # print "Number of bootstraps: ", self.bootstrap
                for i in range(self.bootstrap):
                    # -1 as it is inclusive
                    rdv = randint(0, nr_variants - 1)
                    # print "The random number is: ", rdv,self.variants_sequence[residuenr], self.variants_sequence[residuenr][aa][0], nr_variants, rdv

                    tmp_seq = self.variants_sequence[residuenr][rdv]

                    rdv_val_seq = randint(0, len(self.seq_data[tmp_seq]))
                    tmp_value = float(self.seq_data[tmp_seq][rdv_val_seq])
                    self.bootstrap_results[newkey].append(tmp_value)
            else:

                newkey = residuenr
                if (newkey not in self.bootstrap_results.keys()):
                    self.bootstrap_results[newkey] = []
                rdv = 0
                tmp_seq = self.variants_sequence[residuenr][rdv]

                # print len(self.seq_data[tmp_seq])
                rdv_val_seq = randint(0, len(self.seq_data[tmp_seq]))
                tmp_value = float(self.seq_data[tmp_seq][rdv_val_seq])
                self.bootstrap_results[newkey].append(tmp_value)


        # loop data
        # with open(newkey.replace("_", "") + "_" + self.chain_name + "_" + self.assay + ".csv", 'w') as f:
        print "Debug: ", self.chain_name
        with open(self.chain_name+"data.csv", 'w') as f:
            f.write("Mutation,mean,median,max,min\n")
            for newkey in self.bootstrap_results.keys():
            # print newkey, len(self.bootstrap_results[newkey])
                f.write( newkey+","+str(round( mean(self.bootstrap_results[newkey]),self.round_nr))+","+str(round( median(self.bootstrap_results[newkey]),self.round_nr ))+","+str(round(max(self.bootstrap_results[newkey]),self.round_nr))+","+str(round(min(self.bootstrap_results[newkey]),self.round_nr))+"\n" )


    def max_value(self, inputlist):
        return max([sublist[-1] for sublist in inputlist])

    def min_value(self, inputlist):
        return min([sublist[-1] for sublist in inputlist])

    def get_boxplot(self):
        data = []
        data_names = []
        index = []
        for key in self.bootstrap_results.keys():
            data.append(self.bootstrap_results[key])
            data_names.append( key )
            index.append(int(key.split("_")[0]))


        list1, list2, list3 = (list(x) for x in zip(*sorted(zip(index,data_names, data))))
        data_names = list2
        data = list3

        top = self.max_value(data)
        bottom = self.min_value(data)

        #print "Maximum value: ", top, data[0]
        #print "Minimum value: ", bottom

        fig, ax1 = plt.subplots(figsize=(10,6))
        fig.canvas.set_window_title(self.title)
        #plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
        plt.subplots_adjust(left=0.090, right=0.95 )#, top=0.95, bottom=0.25)

        bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
        plt.setp(bp['boxes'], color='black')
        plt.setp(bp['whiskers'], color='black')
        plt.setp(bp['fliers'], color='red', marker='+')

        # Add a horizontal grid to the plot, but make it very light in color
        # so we can use it for reading data values but not be distracting
        ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5)

        # Hide these grid behind plot objects
        ax1.set_axisbelow(True)
        ax1.set_title(self.title)
        ax1.set_xlabel('Amino Acid Substitution')
        ax1.set_ylabel( self.ylabel )

        # Now fill the boxes with desired colors
        boxColors = ['darkkhaki','royalblue']
        # number of boxes in the plto
        numBoxes = len( data )
        medians = range(numBoxes)


        for i in range(numBoxes):
            box = bp['boxes'][i]
            boxX = []
            boxY = []

            for j in range( 5 ):
                boxX.append(box.get_xdata()[j])
                boxY.append(box.get_ydata()[j])


            boxCoords = zip(boxX,boxY)
            # Alternate between Dark Khaki and Royal Blue
            k = i % 2
            boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
            ax1.add_patch(boxPolygon)

            # Now draw the median lines back over what we just filled in
            med = bp['medians'][i]
            medianX = []
            medianY = []
            for j in range(2):
                medianX.append(med.get_xdata()[j])
                medianY.append(med.get_ydata()[j])
                plt.plot(medianX, medianY, 'k',marker='o')
                medians[i] = medianY[0]
            # Finally, overplot the sample averages, with horizontal alignment
            # in the center of each box
            plt.plot([np.average(med.get_xdata())], [np.average(data[i])],color='w', marker='*', markeredgecolor='k')

        # Set the axes ranges and axes labels
        ax1.set_xlim(0.5, numBoxes+0.5)
        # setting y-limits
        #####ax1.set_ylim(bottom / 1.05, top * 1.10)

        xtickNames = plt.setp(ax1, xticklabels=np.repeat( data_names , 1))
        plt.setp(xtickNames, rotation=45, fontsize=8)

        # Due to the Y-axis scale being different across samples, it can be
        # hard to compare differences in medians across the samples. Add upper
        # X-axis tick labels with the sample medians to aid in comparison
        # (just use two decimal places of precision)
        pos = np.arange(numBoxes)+1
        upperLabels = [str(np.round(s, 2)) for s in medians]
        weights = ['bold', 'semibold']
        for tick,label in zip(range(numBoxes),ax1.get_xticklabels()):
            k = tick % 2
            ##ax1.text(pos[tick], top-(top*0.05), upperLabels[tick],horizontalalignment='center', size='x-small', weight=weights[k],color=boxColors[k])

        plt.axhline(y=self.mean_reference)
        plt.savefig("boxplot_1.png",dpi=600)
        plt.close()

    def main(self):

        parser = argparse.ArgumentParser(description="Bootstrap data to get confidence intervals")
        parser.add_argument('-d', dest='datafile', help='File with sequence and data (Format csv)')
        parser.add_argument('-b', dest='bootstrap', help='Number of bootstraps to be performed',type=int)
        parser.add_argument('-n', dest='nativesequencefile', help='Fasta file with native sequence')
        parser.add_argument('-c', dest='chain_name', help='which chain is being analyzed')
        parser.add_argument('-t', dest='title', help='Title of plot')


        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        # the reference sequence is set
        self.set_native_sequence()

        # set variant data seq to observable
        self.set_variants_data()

        # get the difference between the variants and
        # reference sequence
        self.set_difference()

        # bootstrap values
        self.bootstrap_values()
        #for boot in self.bootstrap_results:
        #    print boot,mean(self.bootstrap_results[boot] ),sqrt(var(self.bootstrap_results[boot] )),self.bootstrap_results[boot]
        self.get_boxplot()

if __name__ == "__main__":
    run = BootstrapData()
    run.main()
