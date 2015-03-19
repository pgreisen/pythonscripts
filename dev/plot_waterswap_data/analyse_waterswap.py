#!/usr/bin/env python
'''
You would need to run two bash-commands to get the data for the analysis. 

see below for the scripts

'''
import sys
from pylab import *
import pylab as plt
import numpy as np
import argparse
from collections import defaultdict
import pdb;

class AnalyseWaterSwap:


    def __init__(self):
        self.burnin = 100
        self.aa_binding_energy = defaultdict(list)
        self.binding_data = ""
        self.aa_binding_data = ""
        self.kd = []
        self.round = 2
        self.boxplot_cutoff = 10
        self.plot_number = 1

    def get_binding_data(self):
        with open(self.binding_data, 'r') as f:
            for line in f:
                tmp = line.split()
                if( int(tmp[5]) >= self.burnin ):
                    self.kd.append( float( line.split()[7] ) )



    def get_aa_binding_data(self):
        with open( self.aa_binding_data, 'r') as f:
            for line in f:
                tmp = line.split()
                key = tmp[1]+"_"+tmp[3]
                self.aa_binding_energy[key].append( float(tmp[5]) )



    def plot_histogram(self):
        mn = round( mean(self.kd), self.round )
        sd = round( sqrt( var(self.kd) ), self.round )
        hist(self.kd)
        xmin,xmax = xlim()
        ymin,ymax = ylim()
        text( xmin + (xmax-xmin)*0.02,  (ymax-ymin)/2,"Mean = "+str(mn)+"\n SD= "+str(sd),size='x-large')
        ## title( filename)
        savefig("binding_constant.png")
        show()
        
    def get_boxplot(self):
        dummy = 0
        data = []
        data_names = []
        for key in self.aa_binding_energy:
            if(key[0:3] != 'WAT'):
                # dummy = 10 make plot
                dummy = dummy + 1
                # print dummy == self.boxplot_cutoff, dummy

                if(dummy <= self.boxplot_cutoff):
                    data.append(self.aa_binding_energy[key])
                    data_names.append(key)

                if( dummy == self.boxplot_cutoff ):

                    fig, ax1 = plt.subplots(figsize=(10,6))
                    fig.canvas.set_window_title('A Boxplot Example')
                    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

                    bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
                    plt.setp(bp['boxes'], color='black')
                    plt.setp(bp['whiskers'], color='black')
                    plt.setp(bp['fliers'], color='red', marker='+')

                    # Add a horizontal grid to the plot, but make it very light in color
                    # so we can use it for reading data values but not be distracting
                    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5)

                    # Hide these grid behind plot objects
                    ax1.set_axisbelow(True)
                    ax1.set_title('Energy Contribution to Binding per Residue')
                    ax1.set_xlabel('Amino Acid')
                    ax1.set_ylabel('Free Energy')
                    # Now fill the boxes with desired colors
                    boxColors = ['darkkhaki','royalblue']
                    # number of boxes in the plto
                    numBoxes = len( data )
                    medians = range(numBoxes)


                    for i in range(numBoxes):
                        box = bp['boxes'][i]
                        boxX = []
                        boxY = []

                        # 11-03-2015
                        # is this the same for all boxes?
                        for j in range( 5 ):

                            # pdb.set_trace()
                            # what does this contain
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
                            plt.plot(medianX, medianY, 'k')
                            medians[i] = medianY[0]
                    # Finally, overplot the sample averages, with horizontal alignment
                    # in the center of each box
                    plt.plot([np.average(med.get_xdata())], [np.average(data[i])],color='w', marker='*', markeredgecolor='k')

                    # Set the axes ranges and axes labels
                    ax1.set_xlim(0.5, numBoxes+0.5)
                    # setting y-limits
                    top = 20
                    bottom = -20
                    ax1.set_ylim(bottom, top)
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
                        ax1.text(pos[tick], top-(top*0.05), upperLabels[tick],horizontalalignment='center', size='x-small', weight=weights[k],color=boxColors[k])

                    # Finally, add a basic legend
                    ##plt.figtext(0.80, 0.08, ' Random Numbers' ,backgroundcolor=boxColors[0], color='black', weight='roman',size='x-small')
                    ##plt.figtext(0.80, 0.045, 'IID Bootstrap Resample',backgroundcolor=boxColors[1],color='white', weight='roman', size='x-small')
                    ##plt.figtext(0.80, 0.015, '*', color='white', backgroundcolor='silver',weight='roman', size='medium')
                    ##plt.figtext(0.815, 0.013, ' Average Value', color='black', weight='roman',size='x-small')
                    ## plt.show()
                    plt.savefig("boxplot_"+str(self.plot_number)+".png")
                    plt.close()
                    self.plot_number = self.plot_number + 1

                    # reset all the plot values:
                    data = []
                    data_names = []
                    dummy = 0



    def main(self):
        parser = argparse.ArgumentParser(description="Analyse results from WaterSwap")
        parser.add_argument('--burnin',dest='burnin', help='Burnin number to skip from the computationa', default=100,type=int )
        parser.add_argument('--binding_data',dest='binding_data', help='Data file with binding data' )
        parser.add_argument('--aa_binding_data',dest='aa_binding_data', help='Data file with aa-binding data' )


        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        self.get_aa_binding_data()
        self.get_boxplot()

        self.get_binding_data()
        self.plot_histogram()


if __name__ == "__main__":
    run = AnalyseWaterSwap()
    run.main()


'''
kd="^Binding";
datapath=/work/greisen/Projects/LigandBinding/Apixaban/WaterSwap_Fentanyl_Apixaban/09032015_Results/binding_costant.dat;
for i in *[0-9]*;
do
    if [ -d "$i" ]; then
        cd $i;
        echo $i;
        # go to the output directory
        cd output;
        for j in results*log;
        do
            grep $kd $j >> $datapath;
        done
        cd ../..;
    fi
done
'''
