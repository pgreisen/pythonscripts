from numpy import mean,sqrt,var
import sys
from collections import defaultdict
from collections import OrderedDict
import os,shutil,sys,argparse

class Analyse2D:

    def __init__(self):
        self.scores = defaultdict(list)
        self.scores_position_two = {}



        self.aa = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        self.aa_single_letter = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

        self.aa_values = []
        # times above mean
        self.factor = 1

        self.score_term = "total_score"
        self.score_term_position = 0
        # A_D_resfile_scores
        self.baseline_value = -537.691


    def get_sorted_hashtable(self, hashtable):
        return OrderedDict(sorted(hashtable.items(), key=lambda x: x[1],reverse=True)) #[0:maxvalue])

    def set_matrix(self,filename):
        tmp_variable = True
        tmp_variable2 = False
        with open(filename,'r') as f:
            for line in f:
                tmp_line = line.split()
                if( line[0:4] == "SEQU" ):
                    continue
                elif ( line[0:4] == "SCOR" and tmp_variable == True):
                    for t in range( len(tmp_line) ):
                        if (tmp_line[t] == self.score_term):
                            self.score_term_position = t
                    tmp_variable = False
                elif( tmp_variable == False):

                    aa_tmp = filename.split('_')
                    tmp_value = round(float( tmp_line[self.score_term_position] ) - self.baseline_value,3)
                    self.scores[ str(aa_tmp[0]) ].append( (aa_tmp[1], str( tmp_value )  ) )

                    # works
                    #self.scores[ str(aa_tmp[0]) ].append( (aa_tmp[1], tmp_line[self.score_term_position] - self.baseline_value)  )

                    ##print  aa_tmp[1] ,  tmp_line[self.score_term_position]
                    ##self.scores[ str(aa_tmp[0]) ].append( scores_position_two[ aa_tmp[1] ] = tmp_line[self.score_term_position]  )

                    #self.scores_position_two = { aa_tmp[1] : tmp_line[self.score_term_position] }
                    #tmp_variable2 = True


                    #elif( tmp_variable2 == True):
                    #self.scores[ str(aa_tmp[0]) ].append( self.scores_position_two )


    def write_matrix_to_file(self):

        with open("2Dscan.csv",'w') as f:
            f.write("AA(203/233),")
            for key in self.aa_single_letter:
                f.write(key+",")

            f.write("\n")
            #assert 1 == 0
            # print self.scores

            for key in self.aa_single_letter:

                f.write(key+",")
                ##print key
                ##assert 1 == 0
                tmp_dic = self.scores[ key ]

                for key in self.aa_single_letter:

                    for i in range(len(tmp_dic) ):

                        if(key == tmp_dic[i][0]):
                            # continue
                            # print key, tmp_dic[i]
                            f.write(tmp_dic[i][1]+",")
                f.write("\n")


    def main(self):

        parser = argparse.ArgumentParser(description="Generate 2D matrix from multiple rosetta output files")

        # get the initial rosetta design as input
        parser.add_argument("-s","--score_term", dest="score_term", help="Which score term to analyse (Default total_score )", default="total_score")

        # parser.add_argument("-b", "--bundles", dest="helical_bundle", help="Four chains helical bundle with four chains is set to true", action="store_true", default=False )

        path = "./"

        files = os.listdir( path )

        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        for fl in files:
            # assumes that the pachdock file ends with .out
            if( os.path.isfile(fl) and fl.endswith("scores") ):
                self.set_matrix( fl )


        self.write_matrix_to_file()


if __name__ == "__main__":
   run = Analyse2D()
   run.main()
