import sys, shutil, os, subprocess, argparse
from pylab import *
from collections import defaultdict
from math import log


class SSMAnalysis:

        def __init__(self):
            self.offset = 0
            self.matrix = defaultdict(list)
            self.datafile = ""
            self.matrix_for_plotting = []
            self.shape_aa = 21
            self.shape_protein = 0
            self.count_sort = []
            self.count_naive = []
            self.column_insert = ["A","V","L","I","M","C","F","Y","W","H","S","T","N","Q","D","E","K","R","G","P","*"]
            # dimensions of figure
            self.x = 12
            self.y = 12
            self.cutoff_count = 10
            self.native_entity = {}



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

                    Asort = float( tmp_line[3] )
                    Bsort = float( tmp_line[4] )

                    # The frequency ration is computed
                    freq = round( log( Asort / Bsort ) ,3 )
                    # freq = round( float(tmp_line[3]) / float(tmp_line[4]) ,3 )

                    # generate a key that contains the native entity as well
                    tmp_resid, resname = tmp_line[0].split('-')
                    key = self.native_entity[tmp_resid]+str(int(tmp_resid)+self.offset)+resname


                    if( Bsort > self.cutoff_count ):
                        self.matrix[ key ] = [ Bsort, Asort , freq ]

                    # store the frequency ratio
                    self.matrix_for_plotting.append( freq )

                    # store the counts from sort
                    self.count_sort.append( Asort )

                    # store the counts from sort
                    self.count_naive.append( Bsort )

            # numpy.asarray
            t = asarray(self.matrix_for_plotting)
            self.shape_protein = len( t ) / self.shape_aa
            self.matrix_for_plotting =  t.reshape( self.shape_protein, self.shape_aa  ).transpose()

            print "Number of residues( Row )  is equal to :", len(self.matrix_for_plotting[:,1])
            print "Amino acids ( column ) is equal to: ", len(self.matrix_for_plotting[1,:])


        def get_native_entity_and_position(self):
            with open( self.datafile,'r') as f:
                for line in f:
                    # 0: is key of residue
                    # 3: A is the sort
                    # 4: B is the naive sort
                    tmp_line = line.split()
                    if( tmp_line[0] == "SeqID"):
                        continue
                    elif( tmp_line[8] == "WT"):
                        key,value = tmp_line[0].split("-")
                        self.native_entity[key] = value

        def plot_matrix(self):
            "Heat map of the matrix"

            #fig = figure(figsize=(self.x,self.y),facecolor='w')
            #ax = fig.add_subplot(111)
            # imshow( self.matrix_for_plotting , interpolation='nearest',cmap='RdYlGn',aspect="equal")
            i = imshow( self.matrix_for_plotting , interpolation='nearest',cmap='RdYlGn', aspect="equal")
            # does not work
            # ax.grid(True)
            # extent=[0,100,0,1], aspect=100

            # fig.colorbar(i )
            # colorbar()
            savefig("heatmap.eps")
            # show()

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

        def get_sorted_frequencies(self):
            """ Return a sorted dictionary based in the set value in this case the frequency
            """
            return sorted( self.matrix.iteritems(),key=lambda (k,v): v[2],reverse=True)


        def main(self):
            parser = argparse.ArgumentParser(description="Analyse SSM data")
            parser.add_argument('-f',dest='datafile', help='File with data from the analysis' )
            parser.add_argument('--offset',dest='offset', default=0,type=int, help="If your sequence is offset by a number one can add or subtract this number")
            parser.add_argument('--cutoff_count',dest='cutoff_count', default=10,type=int, help="Remove enrichments values which have a count less than this in the naive population")
            args_dict = vars( parser.parse_args() )
            for item in args_dict:
                setattr(self, item, args_dict[item])


            self.get_native_entity_and_position()
            self.set_data( )
            self.plot_matrix()
            self.write_matrix_to_csv_file()

            sorted_matrix = self.get_sorted_frequencies()
            with open("frequencies_sorted.txt",'w') as f :
                for key in sorted_matrix:
                    f.write(str(key[0])+","+str(key[1][0])+","+str(key[1][1])+","+str(key[1][2])+"\n" )
                    # f.write(" "+str(key) +"\n"  )


if __name__ == '__main__':
    run = SSMAnalysis()
    run.main()