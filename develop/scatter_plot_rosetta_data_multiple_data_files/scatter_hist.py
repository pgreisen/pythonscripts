import argparse, subprocess,os,sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
'''

The program generates a scatter plot between two data sets.

python ../scatter_hist.py -x sc.dat -y number_of_mutations.dat --xlabel SC --ylabel Mutation


'''

class ScatterHistogram():


    def __init__(self):
        self.xmin = 100
        self.xmax = -1000
        self.ymin = 100
        self.ymax = -1000



    def setup_multiplot(self, xlabel, ylabel ):

        nullfmt   = NullFormatter()         # no labels

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left+width+0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]

        # location of legend
        # does not work
        # loc = [ left_h, bottom_h, 0.2, 0.2 ]


        # start with a rectangular Figure
        plt.figure(1, figsize=(10,10))

        axScatter = plt.axes(rect_scatter)
        axScatter.set_xlabel(xlabel,fontsize=18)
        axScatter.set_ylabel(ylabel,fontsize=18)

        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)

        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)

        return plt, axScatter, axHistx, axHisty

    def scatter_histogram(self, x, y, axScatter, axHistx, axHisty, color_cycle ):

        # the scatter plot:
        axScatter.scatter(x, y, c=color_cycle )

        # now determine nice limits by hand:
        binwidth = 0.25
        xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
        lim = ( int(xymax/binwidth) + 1) * binwidth

        xmin = min(x)
        xmax = max(x)
        ymin = min(y)
        ymax = max(y)

        if( xmin < self.xmin  ):
            self.xmin = xmin
        if( xmax > self.xmax  ):
            self.xmax = xmax
        if( ymin < self.ymin  ):
            self.ymin = ymin
        if( ymax > self.ymax  ):
            self.ymax = ymax

        axScatter.set_xlim( (self.xmin, self.xmax) )
        axScatter.set_ylim( (self.ymin, self.ymax) )

        # bins = np.arange(-lim, lim + binwidth, binwidth)
        # binwidth = (max(x) - min(x))/20
        bins = binwidth*len(x)

        axHistx.hist(x, bins=bins, normed=1, color=color_cycle, stacked=True)
        axHisty.hist(y, bins=bins, orientation='horizontal', normed=1, color=color_cycle, stacked=True)
        
        axHistx.set_xlim( axScatter.get_xlim() )
        axHisty.set_ylim( axScatter.get_ylim() )

        # print self.xmin, self.xmax, self.ymin, self.ymax
        
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
        parser.add_argument("-x", "--xarg", dest="xdata", nargs='*', help="The data has to be paired such that x1 goes with y1 and x2 with y2")
        parser.add_argument("-y", "--yarg", dest="ydata", help="Data set 2", nargs='*')

        # get the data
        parser.add_argument("--xlabel", dest="xlabel", help="Label of x-axis ( First data set) ")
        parser.add_argument("--ylabel", dest="ylabel", help="Label of y-axis ( Second data set")


        parser.add_argument("--label", dest="label", help="Label for each data set", nargs='*', default='None')


        # title of figure
        parser.add_argument("--name", dest="figuretitle", help="Title of figure", default="Title of plot")

        input_variables = parser.parse_args()


        # Number of plots to concatenate
        number_of_plots = len( input_variables.xdata)

        assert number_of_plots == len( input_variables.ydata ), "The number of x data set is different from the y data set"

        # mylegend
        mylegend = []

        color_cycle=['r', 'g', 'b' ] # , 'y'])

        # setup figure
        plt, axScatter, axHistx, axHisty = self.setup_multiplot( input_variables.xlabel, input_variables.ylabel )

        for i in range( number_of_plots ):

            x_data = self.get_data( input_variables.xdata[i] )
            y_data = self.get_data( input_variables.ydata[i] )

            xdata, ydata = self.get_align_vectors( x_data, y_data )

            self.scatter_histogram(xdata ,ydata, axScatter, axHistx, axHisty, color_cycle[i] )

            mylegend.append( input_variables.label[i] )

        # plt.legend( mylegend ) # , loc='upper left')
        plt.savefig( input_variables.figuretitle+".png")
        plt.show()


if __name__ == "__main__":
   run = ScatterHistogram()
   run.main()
