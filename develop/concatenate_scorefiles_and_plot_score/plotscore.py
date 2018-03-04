import pandas as pd
import glob
import argparse,csv
import operator
import seaborn as sns


class ScoreMultiplefiles:


    def __init__(self):
        self.scorex = ""
        self.scorey = ""
        self.scorez = ""
        self.cst_threshold = 7.0
        self.nr_to_file = 5
        self.df = None


    def get_dataframe(self):
        df_ = []
        for scorefile in glob.glob('*.sc'):
            df_.append(pd.read_csv(scorefile,delim_whitespace=True,skiprows=1) )
        return pd.concat(df_)

    def plot_data(self,x,y,dataframe):
        sns_plot = sns.lmplot(x, y, data=dataframe, fit_reg=False)
        sns_plot.savefig("output.png")


    def main(self):

        parser = argparse.ArgumentParser(description="Concatenate scorefiles and plot two value specified and dumb top X number")

        # get the initial rosetta design as input
        parser.add_argument("-x", "--scorex", dest="scorex", help="")
        parser.add_argument("-y", "--scorey", dest="scorey", help="")
        parser.add_argument("-z", "--scorez", dest="scorez", help="")


        args_dict = vars( parser.parse_args() )
        for item in args_dict:
            setattr(self, item, args_dict[item])

        df  = self.get_dataframe()
        stat =  df[self.scorey].describe()
        # print stat.index, type(stat),stat[u'25%'],stat
        print "The 25% value is ", stat[u'25%']

        print df[self.scorez].describe()

        # filter values



        df = df[ ( df[self.scorey] < stat[u'25%']  ) & ( df[self.scorez] < self.cst_threshold  )   ]
        # sort dataframe
        df.sort_values(by=self.scorex,ascending=False, inplace=True)
        self.plot_data(self.scorex,self.scorey,df)
        df.to_csv("lowest_scoring.csv")


if __name__ == "__main__":
   run = ScoreMultiplefiles()
   run.main()




