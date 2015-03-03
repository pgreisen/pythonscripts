from pylab import *
import sys,os
from collections import defaultdict


class get_union:



    def __init__(self):
        self.significant_level = 2
        # self.mutations_positions = {} # defaultdict(list)


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


        mu_data1 = mean( frequencies_scores_dt1 )
        mu_data2 = mean( frequencies_scores_dt2 )

        sd_data1 = sqrt(var(frequencies_scores_dt1))
        sd_data2 = sqrt(var(frequencies_scores_dt2))

        print "Mean of data set 1: ", mu_data1
        print "Standard deviation of set 1: ", sd_data1
        print "Mean of data set 2: ", mu_data2
        print "Standard deviation of set 2: ", sd_data2


        # remove depleted positions
        x = []
        y = []
        labels = []
        for i in intersection:
            # If less than zero skip it:
            if( float(dt1[i]) <= 0.0 and float(dt2[i]) <= 0.0 ):
                continue
            # remove data points which are not significant for both data sets
            elif( float(dt1[i]) <= mu_data1+self.significant_level*sd_data1 and float(dt2[i]) <= mu_data2+self.significant_level*sd_data2 ):
                continue

            else:
                x.append( dt1[i]  )
                y.append( dt2[i]  )
                labels.append(i)
        # f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        # f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)

        scatter(x,y)
        xlim([0, max(frequencies_scores_dt1)])
        ylim([0, max(frequencies_scores_dt2)])

        xlabel("APN")
        ylabel("FEN")

        for label, i, j in zip(labels, x, y):     
            annotate(label, size=12, xy = (i, j), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom', bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.25), arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        savefig("3p_APN_FEN_both_enriched.eps")
        savefig("3p_APN_FEN_both_enriched.png")

        scatter(x,y)
        xlim([0, max(frequencies_scores_dt1)])
        ylim([min(frequencies_scores_dt1), 0])

        xlabel("APN")
        ylabel("FEN")

        for label, i, j in zip(labels, x, y):
            annotate(label, size=12, xy = (i, j), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom', bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.25), arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        savefig("3p_APN_FEN_APNenriched.eps")
        savefig("3p_APN_FEN_APNenriched.png")


        scatter(x,y)
        ylim([1, max(frequencies_scores_dt2)])
        xlim([min( frequencies_scores_dt1), 0])

        xlabel("APN")
        ylabel("FEN")

        for label, i, j in zip(labels, x, y):
            annotate(label, size=12, xy = (i, j), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom', bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.25), arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        savefig("3p_APN_FEN_FENenriched.eps")
        savefig("3p_APN_FEN_FENenriched.png")





if __name__ == "__main__":
    run = get_union()
    run.main()

        

#keys_a = set(dict_a.keys())
#keys_b = set(dict_b.keys())
#intersection = keys_a & keys_b # '&' operator is used for set intersection
