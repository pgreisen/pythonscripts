#from numpy import *
import os,shutil,sys,argparse
import operator
from collections import OrderedDict
from plot_data import *


class analyse_patchdock:

    
    def __init__(self):
        self.scores = {}
        self.score_for_output = {}
        self.debug = 0
        self.path = "./"
        self.cutoff = 2.0
        self.score = 2
        self.mn = 0
        self.sd = 0

    def print_values(self):
        for k, v in self.scores.items():                                                                                    
            print "%s: %s" % (k, v)        
        

    # Requires scorefile 
    # Returns list with lines
    def get_file(self,patchdockfile):
        fl = open(patchdockfile,'r')
        pfile = fl.readlines()
        fl.close()
        return pfile


    # Requires list, maximum_number
    # Returns list with scores in that interval
    def get_scores(self,pfile,filename):
        start = False
        j =  0

        prt = ""
        lig = ""

        for line in pfile:

            tmp_line = line.split()

            length_tmp_line = len( tmp_line )

            if(len(tmp_line) > 0 and  str(tmp_line[0]) == 'Best'):
                break

            # get key for ligand and protein
            if( length_tmp_line == 3 and tmp_line[0] == "ligandPdb" ):
                lig = tmp_line[2][0:3]

                if(self.debug == 1 ):
                    print "The following ligand name is used for the key: ", lig

            # get key for ligand and protein
            if( length_tmp_line == 3 and tmp_line[0] == "receptorPdb" ):
                prt = tmp_line[2][0:4]
                if(self.debug == 1):
                    print "The following ligand name is used for the key: ", prt

            if(start == True):
                j = j + 1
                key = prt+"_"+lig+'_'+str(j)
                # default type is string we have to type cast int
                self.scores[key] = float(tmp_line[ self.score ])

                # get the score out from this - this is temporary fix
                tmp_key = str(filename)+"___"+str(tmp_line[0])
                self.score_for_output [ tmp_key ] = float(tmp_line[ self.score ])


            if(len(tmp_line) > 0 and  str(tmp_line[0]) == '#'):
                start = True

    def get_sorted_results(self): #,maxvalue=10):
        return OrderedDict(sorted(self.scores.items(), key=lambda x: x[1],reverse=True)) #[0:maxvalue])

    def get_sorted_hashtable(self, hashtable): #,maxvalue=10):
        return OrderedDict(sorted(hashtable.items(), key=lambda x: x[1],reverse=True)) #[0:maxvalue])


    def print_number_scores(self):
        print 'The number of scores are ',len(self.scores)

    def get_mean_scores(self):
        return self.mean

    def get_standard_deviation(self):
        return self.sd




    def main(self):
        parser = argparse.ArgumentParser(description="Script to analyse patchdock results")
        parser.add_argument('-f',dest='patchdockscorefiles',nargs='*', help='Specify a number of files to analyse their score')
        parser.add_argument('--path',dest='path', default="./", help='Path to directory where the analysis is performed')
        parser.add_argument('--debug',dest='debug', type=int, help='Run program in debug mode')
        parser.add_argument('--cutoff',dest='cutoff', default=2.0,type=float, help='Display where significant cutoff is')

        parser.add_argument('--score',dest='score', default=2, type=int, help='Score term to analyse - default is 2 which is the total pathdock score')

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        files = os.listdir( self.path )

        for fl in files:
            # assumes that the pachdock file ends with .out
            if( os.path.isfile(fl) and fl.endswith(".out") ):
                tmp_file = self.get_file( fl )
                # we passed the 9 first character as a string to the hash table
                self.get_scores( tmp_file, fl )

        plt = plot_data()
        # get table
        plt.write_table_to_file( self.get_sorted_results() )
        # get mean and standard deviation
        values = plt.get_values( self.get_sorted_results() )

        mn, sd = plt.get_mean( values )

        print "### Mean of data set is: ", round(mn, 4)
        print "### Standard deviation of data set is: ", round(sd, 4)

        if(self.debug == 1):
            print self.cutoff

        plt.plot_histogram_w_mean( values, mn,sd, self.cutoff )

        # get files from dictionary
        cutoff_value = mn+self.cutoff*sd
        print "### The cutoff value for the data set is: ", round(cutoff_value,4)
        plt.write_all_values_above_cutoff( self.get_sorted_hashtable( self.score_for_output ), cutoff_value )


if __name__ == "__main__":
   ana = analyse_patchdock()
   ana.main()



#def main():
#
#    pfile = sys.argv[1]
#    pdbname =  pfile.split('_')[2]
#    max_number = int(sys.argv[2])
#    datafile = get_file(pfile)
#    get_scores(datafile,max_number,pdbname)
#
#    test = sort_results(scores,max_number)
#    
#    for k, v in test.items():
#        print "%s: %s" % (k, v)#
#
#
#if __name__ == "__main__":
#    main()    
################
