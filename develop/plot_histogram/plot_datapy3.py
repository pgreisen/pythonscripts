#!/usr/bin/env python
import os,shutil, sys, math
from numpy import mean
from pylab import *
import argparse,csv
import operator

'''

Assume a single data file first column with names and second column values:


python plot_data_w_mean_sd.py -f FILE --histogram True -n Histogram

'''

class PlotScorefileRosetta:


    def plot_data(self,x,y,filename):
        plot(x,y,'.',c='k',markersize=12)
        # suptitle(name_for_plot)        
        # Added for illustration
        #plot(4.085,-15.733,'o',c='r',markersize=14)
        #plot(9.805,-15.134,'o',c='r',markersize=14)
        savefig(filename+".png")
        show()

    def get_scores_value(self,inputfile, x_string):
        input_file = open(inputfile,'r')
        datafile = []
        x = []
        i = 0
        for line in input_file:
            # Skip the first line
            
            if( line[0:4] == "SEQU" ):
                pass

            elif( line[0:4] == "SCOR" and i != 1 and line.split()[1] == "total_score"):
                
                continue

            elif( line[0:4] == "SCOR" and i == 1):
                x_index = 0
                indexer = 0
                tmpline = line.split()

                for score in tmpline:

                    if ( x_string == score ):

                        x_index = indexer

                    indexer += 1

                assert x_index != 0

            elif( line[0:4] == "SCOR" and i != 1 ):

                tmpline =  line.split()

                x.append( float( tmpline[ x_index ] ) )

            i += 1

        input_file.close()

        return x



    def plot_histogram(self, data, filename):
        hist( data)
        xmin,xmax = xlim()
        ymin,ymax = ylim()
        mn, sd = self.get_mean_and_standard_deviation(data)
        text( xmin + (xmax-xmin)*0.02,  (ymax-ymin)/2,"Mean = "+str(mn)+"\n SD= "+str(sd),size='x-large')
        title( filename)
        savefig("histogram_"+filename+".png")




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

            elif( line[0:4] == "SCOR" and i != 1 and line.split()[1] == "total_score"):
                continue
                
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

                if( len(tmpline) < tag_index ):
                    continue

                x.append( tmpline[ x_index ] )
                y.append( tmpline[ y_index ] )
                tag.append( tmpline[ tag_index ] )

            i += 1

        input_file.close()

        return x,y,tag


    def write_csv_table_rmsd_analysis(self,tagstring, xstring, ystring):

        p_y = []
        p_x = []

        with open('datafile.csv', 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            #spamwriter.writerow(["Tag",str(xstring),str(ystring)])

            length_tag = len( tagstring )
            # keydict = {}
            assert length_tag == len( xstring ) == len( ystring )
            dummy = 1
            for index in range(length_tag):
                spamwriter.writerow([ tagstring[index], xstring[index], ystring[index] ])
                # keydict[ tagstring[index]+" "+str(dummy) ] = float(ystring[index] )
                p_x.append( tagstring[index]+":  "+xstring[index] )
                p_y.append( tagstring[index]+":  "+ystring[index] )

                dummy += 1

        # sorted_x = sorted(keydict.items(), key=operator.itemgetter(1),reverse=True)
        # sorted_x = sorted(keydict.items(), key=operator.itemgetter(1) )
        with open("px.dat",'w') as f:
            for key in p_x:
                # Only get poses where the total energy is less than zero
                f.write( str( key+"\n" ) )
        with open("py.dat",'w') as f:
            for key in p_y:
                # Only get poses where the total energy is less than zero
                f.write( str( key+"\n" ) )


    def get_mean_and_standard_deviation(self,data):
        return round( mean(data), 3 ), round( sqrt( var(data) ), 3)


    def get_data_single_column( self, file ):
        x = []
        input_file = open( file,'r')
        for line in input_file:
            x.append( float( line.split()[0] ) )
        return x


    def main(self):
        print("How to run the script")
        print("Argument 1: Rosetta score file")
        print("Argument 2: Data to collect for the x-axis ( Default=ligand_rms_no_super_B)")
        print("Argument 3: Data to collect for the y-axis ( Default=interface_delta_B)")
        print("Argument 4: Name of plot.png ( Default=results)")


        parser = argparse.ArgumentParser(description="Plot and analyse data from Rosetta score file. Does a 2-D plot ( x versus y) ")

        # get the initial rosetta design as input
        parser.add_argument("-f", "--file", dest="scorefile", help="Initial Rosetta score file")

        parser.add_argument("-x", "--xstring", dest="x_string", help="Score term for the X-axis", default="ligand_rms_no_super_X")

        parser.add_argument("-y", "--ystring", dest="y_string", help="Score term for the Y-axis", default="interface_delta_X")

        parser.add_argument("-n", "--name", dest="filename", help="Name for the plot", default="results")

        parser.add_argument("--histogram", dest="histogram", help="Only plot one varaible", default=0, type=int )

        parser.add_argument("-c", dest="column", help="Only plot one varaible", default="False", type=bool)

        input_variables = parser.parse_args()


        if( input_variables.histogram == 0 ):
            x,y,tag = self.get_scores(input_variables.scorefile,input_variables.x_string, input_variables.y_string)

            self.plot_data(x,y, input_variables.filename)
            self.write_csv_table_rmsd_analysis(tag, x, y)


        elif( input_variables.histogram == 1 and input_variables.column == True):
            data = self.get_data_single_column( input_variables.scorefile )
            self.plot_histogram( data, input_variables.filename )
            self.get_mean_and_standard_deviation( data )


        else:
            data = self.get_scores_value(input_variables.scorefile, input_variables.x_string )
            self.plot_histogram( data, input_variables.filename )
            self.get_mean_and_standard_deviation( data )


if __name__ == "__main__":
   run = PlotScorefileRosetta()
   run.main()
