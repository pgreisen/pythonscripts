#__author__="pgreisen"
import argparse, subprocess,os,sys
from numpy import *
# importing the class
from GetDataScoreFiles import *

class IterativeDesignDocking:

    def __init__(self):
        self.cluster_tag = [] #array([])
        self.cluster_rmsd = [] #array([])
        self.cluster_ife = [] #array([])

    def get_data(self,filename):
        tag = [] #array([])
        rmsd = [] #array([])
        ife = [] # array([])
        with open(filename,'r') as f:
            for line in f:
                tmpline = line.split(',')
                tag.append( tmpline[0])
                rmsd.append( float( tmpline[1]) )
                ife.append( float ( tmpline[2]) )
        return tag, rmsd, ife

    def get_lowest_scoring_point(self,energy_array):
        return argmin( energy_array )

    def cluster_results(self,tag, rmsd, ife, reu_cutoff, rmsd_cutoff):
        number_of_docks = len(tag)
        assert number_of_docks == len(rmsd) == len(ife)

        lowest_scoring_value_position = self.get_lowest_scoring_point( ife )

        # Only collect structure greater than the rmsd cutoff

        if( rmsd[ lowest_scoring_value_position] >= rmsd_cutoff ):

            self.cluster_tag.append(tag[ lowest_scoring_value_position])
            self.cluster_ife.append(ife[ lowest_scoring_value_position])
            self.cluster_rmsd.append(rmsd[ lowest_scoring_value_position])
        iterative = True
        ife_cutoff = ife [lowest_scoring_value_position ] + reu_cutoff
        iter_rmsd_cutoff = rmsd_cutoff


        tag.pop( lowest_scoring_value_position )
        ife.pop( lowest_scoring_value_position)
        rmsd.pop(lowest_scoring_value_position)


        while (iterative):
            lowest_scoring_value_position = self.get_lowest_scoring_point( ife )

            if( ife[ lowest_scoring_value_position  ] > ife_cutoff ):
                iterative = False
                break
            if ( rmsd[ lowest_scoring_value_position ] - iter_rmsd_cutoff > 0.0 ):

                iter_rmsd_cutoff = iter_rmsd_cutoff + rmsd[ lowest_scoring_value_position ]
                self.cluster_tag.append(tag[ lowest_scoring_value_position ])
                self.cluster_ife.append(ife[ lowest_scoring_value_position ])
                self.cluster_rmsd.append(rmsd[ lowest_scoring_value_position ])


                tag.pop( lowest_scoring_value_position )
                ife.pop( lowest_scoring_value_position)
                rmsd.pop(lowest_scoring_value_position)
            else:

                tag.pop( lowest_scoring_value_position )
                ife.pop( lowest_scoring_value_position)
                rmsd.pop(lowest_scoring_value_position)

                continue

        return self.cluster_tag, self.cluster_ife, self.cluster_rmsd

    # setup design runs
    def setup_design(self):
        pass


    def main(self):

        parser = argparse.ArgumentParser(description='IterativeDesignDocking - starting from a docking run')

        parser.add_argument("-f", "--file", dest="pdbfile", help="File with aligned ligand")

        parser.add_argument("-s", "--scorefile", dest="scorefile", help="CSV file with tag, RMSD, IFE")

        parser.add_argument( "--reu", dest="reu_cutoff", help="REU within cutoff - allow e.g. 2 REU clusters", type=float,
                             default=1.0)

        parser.add_argument( "--rmsd", dest="rmsd_cutoff", help="RMSD within cutoff - allow e.g. 2 RMSD clusters", type=float,
                             default=1.0)

        parser.add_argument( "-x", dest="x_string", help="String containing the score type you want for one dimension (Default=ligand_auto_rms_no_super)", type=str,
                             default="ligand_auto_rms_no_super")

        parser.add_argument( "-y", dest="y_string", help="String containing the score type you want for one dimension (Default=interface_delta)", type=str,
                             default="interface_delta")

        parser.add_argument( "--scoretype", dest="scoretype", help="What kind of score file is used as input ( default=score.sc )", type=str,
                             default="sc")

        input_arguments = parser.parse_args()

        if( input_arguments.scoretype == "sc"):
            gdsf = GetDataScoreFiles()
            # return x_string, y_string, tags
            rmsd, ife, tag = gdsf.get_scores( input_arguments.scorefile, input_arguments.x_string, input_arguments.y_string )
        elif( input_arguments.scoretype == "csv" ):
            tag, rmsd, ife = self.get_data( input_arguments.scorefile )
        else:
            print "Unknown scoretype is given"
            sys.exit()


        # lowest_scoring_position_in_array
        cluster_tag, cluster_rmsd, cluster_ife = self.cluster_results( tag, rmsd, ife, input_arguments.reu_cutoff, input_arguments.rmsd_cutoff)


        if( len( cluster_tag ) == 0 ):
            print "The input structure has the lowest energy and rmsd hence no further improvements are needed"
            sys.exit("%%%%%%%%%%%%%% Script is ending %%%%%%%%%%%%%%%%%%\n")

        print "The following poses score lower than the native",cluster_tag

if __name__ == "__main__":
   run = IterativeDesignDocking()
   run.main()