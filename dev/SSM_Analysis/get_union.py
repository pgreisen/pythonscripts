from pylab import *
import sys,os
from matplotlib import rc
from collections import defaultdict


class get_union:



    def __init__(self):
        self.significant_level = 2
        self.round = "round4"
        self.prime = "3p"
        # self.mutations_positions = {} # defaultdict(list)
        self.cutoff_min_x = float( log(10) )



    def get_data(self,datafile):
        hashtable = {}

        frequencies_score = []

        with open(datafile,'r') as f:
            for line in f:
                tmp = line.split(',')
                hashtable[tmp[0]] = tmp[3]
                frequencies_score.append(float(tmp[3]))
        return hashtable, frequencies_score


    def get_intersection_values(self,intersection,dt1,dt2):
        values = {}
        for i in intersection:
            values[ dt1[i] ] = i+"Rosetta"
            values[ dt2[i] ] = i+"FoldX"
        return values

    def main(self):

        datafile_1 = sys.argv[1]
        datafile_2 = sys.argv[2]

        dt1, frequencies_scores_dt1 = self.get_data( datafile_1 )
        dt2, frequencies_scores_dt2 = self.get_data( datafile_2 )

        keys_a = set(dt1.keys())
        keys_b = set(dt2.keys())
        intersection = keys_a & keys_b # '&' operator is used for set intersection

        # unique_keys_a =  list(set( keys_a ) - set(keys_b) )

        mu_data1 = mean( frequencies_scores_dt1 )
        mu_data2 = mean( frequencies_scores_dt2 )

        sd_data1 = sqrt(var(frequencies_scores_dt1))
        sd_data2 = sqrt(var(frequencies_scores_dt2))

        print "Mean of data set 1: ", mu_data1
        print "Standard deviation of set 1: ", sd_data1
        print "Mean of data set 2: ", mu_data2
        print "Standard deviation of set 2: ", sd_data2


        # remove depleted positions
        # positions are enriched for both positions
        x1 = []
        y1 = []
        labels1 = []

        # positions are enriched for Fentanyl but not for Apixaban
        x2 = []
        y2 = []
        labels2 = []

        # positions are enriched for Apixaban but not for Fentanyl
        x3 = []
        y3 = []
        labels3 = []



        for i in intersection:

            # If less than zero skip it:
            if( float(dt1[i]) <= 0.0 and float(dt2[i]) <= 0.0 ):
                continue

            # remove data points which are not significant for both data sets
            elif( float(dt1[i]) <= mu_data1+self.significant_level*sd_data1 and float(dt2[i]) <= mu_data2+self.significant_level*sd_data2 ):
                continue

            # both positions are enriched
            elif( float(dt1[i]) >= self.cutoff_min_x and float(dt2[i]) >= self.cutoff_min_x ):
                x1.append( dt1[i] )
                y1.append( dt2[i] )
                labels1.append( i )

           # both positions are enriched
            elif( float(dt1[i]) >= self.cutoff_min_x and float(dt2[i]) < self.cutoff_min_x):

                print dt1[i], dt2[i]
                x2.append( dt1[i] )
                y2.append( dt2[i] )
                labels2.append( i )

          # both positions are enriched
            elif( float(dt2[i]) >= self.cutoff_min_x and float(dt1[i]) < self.cutoff_min_x):
                x3.append( dt1[i] )
                y3.append( dt2[i] )
                labels3.append( i )

            else:
                # print i, dt1[i],dt2[i]
                continue
                #


        # f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        # f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        # rc('text', usetex=True)



        ##title("3' Sort 4 Enrichment Fentanyl/Apixaban",size=22)
        scatter(x1,y1,color='red',s=40)
        #xlim([2.6, 3.2])
        ylim([4, 5])
        xlabel(r"Apixaban log$_{2}$(Sort4/Naive)",size=16)
        ylabel(r"Fentanyl log$_{2}$(Sort4/Naive)",size=16)

        for label, i, j in zip(labels1, x1, y1):
            annotate(label, size=12, xy = (i, j), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom', bbox = dict(boxstyle = 'round,pad=0.5', fc = 'gray', alpha = 0.25), arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        savefig(self.prime+"_APN_FEN_both_enriched"+self.round+".eps")
        savefig(self.prime+"_APN_FEN_both_enriched"+self.round+".png")
        close()

        ##title("3' Sort 4 Enrichment Apixaban")
        scatter(x2,y2, color='red', s=40)
        #xlim([self.cutoff_min_x , max(frequencies_scores_dt1)])
        #ylim([min(frequencies_scores_dt1), 0])

        xlabel(r"Apixaban log$_{2}$(Sort4/Naive)",size=16)
        ylabel(r"Fentanyl log$_{2}$(Sort4/Naive)",size=16)

        for label, i, j in zip(labels2, x2, y2):
            annotate(label, size=12, xy = (i, j), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom', bbox = dict(boxstyle = 'round,pad=0.5', fc = 'gray', alpha = 0.25), arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        savefig(self.prime+"_APN_FEN_APNenriched"+self.round+".eps")
        savefig(self.prime+"_APN_FEN_APNenriched"+self.round+".png")
        close()

        ##title("3' Sort 4 Enrichment Fentanyl")
        scatter(x3,y3,color='red',s=40)
        #ylim([self.cutoff_min_x , max(frequencies_scores_dt2)])
        #xlim([min( frequencies_scores_dt1), 0])

        xlabel(r"Apixaban log$_{2}$(Sort4/Naive)",size=16)
        ylabel(r"Fentanyl log$_{2}$(Sort4/Naive)",size=16)

        for label, i, j in zip(labels3, x3, y3):
            annotate(label, size=12, xy = (i, j), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom', bbox = dict(boxstyle = 'round,pad=0.5', fc = 'gray', alpha = 0.25), arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        savefig(self.prime+"_APN_FEN_FENenriched"+self.round+".eps")
        savefig(self.prime+"_APN_FEN_FENenriched"+self.round+".png")
        close()

if __name__ == "__main__":
    run = get_union()
    run.main()
