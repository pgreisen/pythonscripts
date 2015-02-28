import sys, shutil, os, subprocess, argparse
from pylab import *
from collections import defaultdict
from math import log


class SSMAnalysis:

        def __init__(self):
            self.offset = 0
            self.matrix = defaultdict(list)
            self.datafile = ""
            self.matrix_for_plotting = None
            self.matrix_for_plotting = []
            self.shape_aa = 21
            self.shape_protein = 0
            self.count_sort = []
            self.count_naive = []
            self.column_insert = ["A","V","L","I","M","C","F","Y","W","H","S","T","N","Q","D","E","K","R","G","P","*"]
            # dimensions of figure
            self.x = 12
            self.y = 12


        def set_data(self ):
            """Reads the data from the input file and sets the values in the hash-table"""
            with open( self.datafile,'r') as f:
                for line in f:
                    # 0: is key of residue
                    # 3: A is the sort
                    # 4: B is the naive sort
                    tmp_line = line.split()
                    if( tmp_line[0] == "SeqID"):
                        continue
                    # The frequency ration is computed
                    freq = round( log( float(tmp_line[3]) / float(tmp_line[4]) ),3 )
                    # freq = round( float(tmp_line[3]) / float(tmp_line[4]) ,3 )
                    self.matrix [tmp_line[0] ].append( ( tmp_line[3], tmp_line[4], freq ) )

                    # store the frequency ratio
                    self.matrix_for_plotting.append( freq )

                    # store the counts from sort
                    self.count_sort.append( float(tmp_line[3]) )

                    # store the counts from sort
                    self.count_naive.append( float(tmp_line[4]) )

            # numpy.asarray
            t = asarray(self.matrix_for_plotting)
            self.shape_protein = len( t ) / self.shape_aa
            self.matrix_for_plotting =  t.reshape( self.shape_protein, self.shape_aa  ).transpose()


            print "ROW IS EQUAL TO :", len(self.matrix_for_plotting[:,1])
            print "Column is equal to: ", len(self.matrix_for_plotting[1,:])



        def plot_matrix(self):
            "Heat map of the matrix"

            fig = figure(figsize=(self.x,self.y),facecolor='w')
            ax = fig.add_subplot(111)
            # imshow( self.matrix_for_plotting , interpolation='nearest',cmap='RdYlGn',aspect="equal")
            i = imshow( self.matrix_for_plotting , interpolation='nearest',cmap='RdYlGn', aspect="equal")
            # does not work
            # ax.grid(True)
            # extent=[0,100,0,1], aspect=100

            # fig.colorbar(i )
            # colorbar()
            show()

        def get_data(self):
            """returns the matrix """
            return self.matrix

        def write_matrix_to_csv_file(self):
            dummy = 0
            with open("raw_rations.dat",'w') as f:
                # write residue number in the first line of the csv
                f.write("")
                for i in range( self.shape_protein ):
                    f.write( str(i+self.offset)+",")
                f.write("\n")
                r_sort = asarray(self.count_sort).reshape( self.shape_aa, self.shape_protein )
                naive_sort = asarray(self.count_naive).reshape( self.shape_aa, self.shape_protein )
                for j in range( self.shape_aa):
                    f.write( str(self.column_insert[dummy])+",")
                    dummy += 1
                    for i in range( self.shape_protein ):
                        f.write(str(r_sort[j][i])+"/"+str(naive_sort[j][i])+"," )
                    f.write("\n")




        def main(self):
            parser = argparse.ArgumentParser(description="Analyse SSM data")
            parser.add_argument('-f',dest='datafile', help='File with data from the analysis' )
            parser.add_argument('--offset',dest='offset', default=0,type=int, help="If your sequence is offset by a number one can add or subtract this number")
            args_dict = vars( parser.parse_args() )
            for item in args_dict:
                setattr(self, item, args_dict[item])

            self.set_data( )
            self.plot_matrix()
            self.write_matrix_to_csv_file()

if __name__ == '__main__':
    run = SSMAnalysis()
    run.main()