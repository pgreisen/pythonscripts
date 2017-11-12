#!/usr/bin/env python
import os,shutil, commands, sys, math
from numpy import mean
from pylab import *
import argparse,csv
from operator import itemgetter
import subprocess

class PlotScorefileRosetta:


    def __init__(self):
        self.number_of_pdbs_to_dump = 50
        self.jobname = ""
        self.datafilecsv = "datafile.csv"
        self.silentout = "silent.out"

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
        i = 1


        for line in input_file:
            # Skip the first line
            if( line[0:4] == "SEQU" ):
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

            elif( line[0:4] == "SCOR" and i != 1 ):

                tmpline =  line.split()

                x.append( tmpline[ x_index ] )
                y.append( tmpline[ y_index ] )
                tag.append( tmpline[ tag_index ] )

            i += 1

        input_file.close()

        return x,y,tag        


    def write_csv_table_rmsd_analysis(self,tagstring, xstring, ystring):

        # write the top scoring docks to file
        top_docks = open("get_top_scored_docks.sh", "w")

        dummy = 1

        # test 26-06-2014
        sorted_by_ystring = sorted(range(len(ystring)), key=lambda x: float(ystring[x]) )
        
        #sorted_by_ystring = sorted( ystring, key=float , reverse=True)
        # a, key=lambda x: float(x))

        with open(self.datafilecsv, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

            length_tag = len( tagstring )

            assert length_tag == len( xstring ) == len( ystring )

            assert len( sorted_by_ystring ) == length_tag
            for index in sorted_by_ystring:
                spamwriter.writerow([ tagstring[index], xstring[index], ystring[index] ])
                if( dummy < self.number_of_pdbs_to_dump ):
                    tmppdb = "~greisen/git_rosetta/rosetta_master/Rosetta/main/source/bin/extract_atomtree_diffs.default.linuxgccrelease -extra_res_fa $1 -in:file:silent "+self.silentout+" -in::file::silent_struct_type binary @/work/greisen/files/rescore_w_facts/flags -database ~greisen/git_rosetta/rosetta_master/Rosetta/main/database -tags "+str(tagstring[index])
                    top_docks.write(tmppdb+"\n")
                    dummy += 1
            top_docks.close()


    def extract_top_dock_poses(self, ):
        '''
        
        '''



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

        parser.add_argument("-x", "--xstring", dest="x_string", help="Score term for the X-axis", default="ligand_rms_no_super_X")

        parser.add_argument("-y", "--ystring", dest="y_string", help="Score term for the Y-axis", default="interface_delta_X")


        parser.add_argument("-n", "--name", dest="filename", help="Name for the plot", default="results")


        parser.add_argument("-j", "--jobname", dest="jobname", help="name of job", default="")


        parser.add_argument("-l", "--ligandparameter", dest="ligandpath", help="Extract 10 lowest scoring models", default="")

        input_variables = parser.parse_args()


        self.jobname = input_variables.jobname

        if (len(self.jobname) > 0):
            self.datafilecsv = self.jobname+"_datafile.csv"
            self.silentout = self.jobname+"_silent.out"

        x,y,tag = self.get_scores(input_variables.scorefile,input_variables.x_string, input_variables.y_string)

        self.plot_data(x,y, input_variables.filename)

        self.write_csv_table_rmsd_analysis(tag, x, y)

        if( len( input_variables.ligandpath ) > 0):
            tmpstring = "sh get_top_scored_docks.sh "+str(input_variables.ligandpath)
            subprocess.Popen( tmpstring,shell=True).wait()


if __name__ == "__main__":
   run = PlotScorefileRosetta()
   run.main()
