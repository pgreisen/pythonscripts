#!/usr/bin/env python
import os,shutil, commands, sys, math
from numpy import mean
from pylab import *
import argparse,csv


class GetDataScoreFiles:


    def plot_data(self,x,y,filename):
        plot(x,y,'.',c='k')
        # suptitle(name_for_plot)
        savefig(filename+".png")

    def get_scores(self,inputfile,x_string, y_string):
        input_file = open(inputfile,'r')
        datafile = []
        x = []
        y = []
        tag = []
        i = 0
        for line in input_file:
            # Skip the first line
            if( line[0:4] == "SEQU" ):
                pass
            elif( line[0:4] == "SCOR" and i == 1):
                indexer = 0
                x_index = 0
                y_index = 0
                tag_index = 0
                tmpline = line.split()
                for score in tmpline:
                    if ( x_string == score ):
                        x_index = indexer
                    elif( y_string == score ):
                        y_index = indexer
                    elif( "description" == score ):
                        tag_index = indexer
                    indexer += 1
                assert tag_index != x_index != y_index
            elif( line[0:4] == "SCOR" and i != 1 ):
                tmpline =  line.split()
                x.append( float( tmpline[ x_index ] ) )
                y.append( float (tmpline[ y_index ] ) )
                tag.append( tmpline[ tag_index ] )
            i += 1
        input_file.close()
        return x,y,tag


    def write_csv_table_rmsd_analysis(self,tagstring, xstring, ystring):
        with open('datafile.csv', 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(["Tag",str(xstring),str(ystring)])

            length_tag = len( tagstring )

            assert length_tag == len( xstring ) == len( ystring )
            for index in range(length_tag):
                spamwriter.writerow([ tagstring[index], xstring[index], ystring[index] ])



    def get_mean_and_standard_deviation(self,data):

        print "The mean of the data set is ", round( mean(data), 3 )
        print "The standard deviation is ", round( sqrt( var(data) ), 3)



    def main(self):

        print "How to run the script"
        print "Argument 1: Rosetta score file"
        print "Argument 2: Data to collect for the x-axis ( Default=ligand_rms_no_super_B)"
        print "Argument 3: Data to collect for the y-axis ( Default=interface_delta_B)"
        print "Argument 4: Name of plot.png ( Default=results)"


        parser = argparse.ArgumentParser(description="Plot and analyse data from Rosetta score file. Does a 2-D plot ( x versus y) ")

        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="scorefile", help="Initial Rosetta score file")

        parser.add_argument("-x", "--xstring", dest="x_string", help="Score term for the X-axis", default="ligand_rms_no_super_B")

        parser.add_argument("-y", "--ystring", dest="y_string", help="Score term for the Y-axis", default="interface_delta_B")


        parser.add_argument("-n", "--name", dest="filename", help="Name for the plot", default="results")

        input_variables = parser.parse_args()


        x,y,tag = self.get_scores(input_variables.scorefile,input_variables.x_string, input_variables.y_string)

        self.plot_data(x,y, input_variables.filename)

        self.write_csv_table_rmsd_analysis(tag, x, y)

        # get_mean_and_standard_deviation(data)


if __name__ == "__main__":
   run = PlotScorefileRosetta()
   run.main()
