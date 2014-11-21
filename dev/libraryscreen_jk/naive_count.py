#!/usr/bin/python
from optparse import OptionParser
import argparse
from pylab import *
from numpy import *

class GetStats():

    def __init__(self):
        self.experimental_data = []


    def collect_data(self,datafile):
        with open(datafile,'r') as f:
            for line in f:
                self.experimental_data.append( float( line.split()[0] ) )

    def plot_data(self):
        bp = boxplot( self.experimental_data)

        print "The median of the data set", bp['medians'][0].get_ydata()
        print "The mean of the data set ", mean( self.experimental_data)
        print "The standard deviation of the data set ", sqrt(var( self.experimental_data))
        print "Maximum value of data set", max( self.experimental_data )


        show()



    def main(self):
        parser = OptionParser()
        parser = argparse.ArgumentParser(description='Analysis a raw count file')
        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="datafile", help="datafile")
        inputs = parser.parse_args()
        self.collect_data( inputs.datafile )

        self.plot_data()


if __name__ == '__main__':
    run = GetStats()
    run.main()

