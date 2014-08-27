import argparse, subprocess,os,sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
'''

The program generates a scatter plot between two data sets.

python ../scatter_hist.py -x sc.dat -y number_of_mutations.dat --xlabel SC --ylabel Mutation


'''

class ScatterHistogram():


    def scatter_histogram(self,x,y,name_of_plot,xlabel,ylabel):

        nullfmt   = NullFormatter()         # no labels
    
        # definitions for the axes 
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left+width+0.02
        
        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]
        
        # start with a rectangular Figure
        plt.figure(1, figsize=(8,8))

        axScatter = plt.axes(rect_scatter)
        axScatter.set_xlabel(xlabel,fontsize=18)
        axScatter.set_ylabel(ylabel,fontsize=18)

        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        
        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        # the scatter plot:
        axScatter.scatter(x, y)

        # now determine nice limits by hand:
        binwidth = 0.25
        xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
        lim = ( int(xymax/binwidth) + 1) * binwidth
        
        xmin = min(x)
        xmax = max(x)
        ymin = min(y)
        ymax = max(y)

        axScatter.set_xlim( (xmin, xmax) )
        axScatter.set_ylim( (ymin, ymax) )

        # bins = np.arange(-lim, lim + binwidth, binwidth)
        # binwidth = (max(x) - min(x))/20
        bins = binwidth*len(x)


        axHistx.hist(x, bins=bins)
        axHisty.hist(y, bins=bins, orientation='horizontal')
        
        axHistx.set_xlim( axScatter.get_xlim() )
        axHisty.set_ylim( axScatter.get_ylim() )
        plt.savefig(name_of_plot+".png")
        plt.show()
        
    def get_data( self, data):
        dict_data = {}
        tmpfile = open(data,'r')
        for line in tmpfile:
            key,value = line.split()
            
            if( ":" in key ):
                key = key.split(':')[0]

            dict_data[key] = value
        return dict_data

    def get_align_vectors(self,a,b):
        """ Will align the two vectors """
        x_values = []
        y_values = []
        for k, v in a.items():
            x_values.append( float( v))
            y_values.append( float( b[k]))

        return x_values, y_values


    def main( self ):

        parser = argparse.ArgumentParser(description='Analyse two data set and plot their correlation')

        # get the data
        parser.add_argument("-x", "--xarg", dest="xdata", help="Data set 1")
        parser.add_argument("-y", "--yarg", dest="ydata", help="Data set 2")

        # get the data
        parser.add_argument("--xlabel", dest="xlabel", help="Label of x-axis ( First data set) ")
        parser.add_argument("--ylabel", dest="ylabel", help="Label of y-axis ( Second data set")

        # title of figure
        parser.add_argument("--name", dest="figuretitle", help="Title of figure", default="Title of plot")

        input_variables = parser.parse_args()

        x_data = self.get_data( input_variables.xdata )
        y_data = self.get_data( input_variables.ydata )

        xdata, ydata = self.get_align_vectors( x_data, y_data )

        self.scatter_histogram(xdata ,ydata , input_variables.figuretitle, input_variables.xlabel, input_variables.ylabel )

if __name__ == "__main__":
   run = ScatterHistogram()
   run.main()
