#!/usr/bin/env python
import os,shutil, commands
import sys, math
from numpy import mean
from pylab import *    

def plot_gnuplot(datafile,title='InsertTitle'):
    #    print datafile,"datafile",type(datafile)
    # make gnuplot input file
    fname = 'plot'

    f = open( fname + '.gnuplot' , 'w')

    f.write("""
    #set title '%s' """ % (title))

    f.write("""
    # Make a png plot
    set term png small;
    set output '%s.png';
    set xlabel 'RMSD to initial position';
    set ylabel 'Interface Energy';
    # basic plot command
    plot '%s' using 1:2 title 'Pose';
    """ % (title,datafile))
    f.close()
    cmd = 'gnuplot -geometry 800x200 -persist '+ fname + '.gnuplot'
    failure,output = commands.getstatusoutput( cmd )
    if failure:
        print 'running gnuplot failed\n%s\n%s' %( cmd, output);
        sys.exit(1)

def plot_histogram( data, filename):
    hist( data)
    mn, sd = get_mean_and_standard_deviation(data)
    xmin,xmax = xlim()
    ymin,ymax = ylim()
    text( xmin + (xmax-xmin)*0.02,  (ymax-ymin)/2,"Mean = "+str(mn)+"\n SD = "+str(sd),size="x-large")
    title( filename)
    savefig("histogram_"+filename+".png")

def read_inputfile(inputfile):
    tmp_file = open(inputfile,'r')
    datafile = []
    for line in tmp_file:
        if( line[0] == '@' or line[0] == '&'):
            continue
        datafile.append( float( line.split()[7]) )
    tmp_file.close()
    return datafile

def get_mean_and_standard_deviation(data):
    mn = round( mean(data), 3 )
    sd = round( sqrt( var(data) ), 3)
    print "The mean of the data set is ", round( mean(data), 3 )
    print "The standard deviation is ", round( sqrt( var(data) ), 3)
    return mn, sd


def main():
    datafile = sys.argv[1]

    data = read_inputfile(datafile)
    
    name_for_plot = sys.argv[2]

    plot_histogram(data,name_for_plot)

    get_mean_and_standard_deviation(data)


if __name__ == "__main__":
    main()


