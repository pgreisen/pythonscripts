#from numpy import *
import os,shutil,sys
import operator
from collections import OrderedDict

class analyse_patchdock:

    
    def __init__(self):
        print "Analysing results from patchdock"
        self.scores = {}

    #    def __str__(self):
    #        print self.scores

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
    def get_scores(self,pfile,pdbname):
        print "Getting score"
        start = False
        j =  0
        for line in pfile:
            tmp_line = line.split()

            if(len(tmp_line) > 0 and  str(tmp_line[0]) == 'Best'):
                break

            if(start == True):
                j = j + 1
                key = pdbname+'_'+str(j)
                # default type is string we have to type cast it
                self.scores[key] = int(tmp_line[2])

            if(len(tmp_line) > 0 and  str(tmp_line[0]) == '#'):
                start = True

    def get_sorted_results(self): #,maxvalue=10):
        return OrderedDict(sorted(self.scores.items(), key=lambda x: x[1],reverse=True)) #[0:maxvalue])


    def print_number_scores(self):
        print 'The number of scores are ',len(self.scores)


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
