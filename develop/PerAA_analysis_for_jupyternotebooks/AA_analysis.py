import pandas as pd
import numpy as np
from scipy.stats.mstats import gmean
import plotly
from plotly.graph_objs import Scatter, Layout
import plotly
import plotly.offline as offline
import plotly.graph_objs as go
import plotly.plotly as py
import plotly.figure_factory as ff
from datetime import date

class AA_analysis:


    def __init__(self):
        self.nr_bootstraps = 750
        # name of column in df that contains mutation name
        self.mutations = "Mutations"
        self.cutoff = 1
        self.today = str(date.today())


    def set_values(self, muts, ratio, df_ratios):
        mt = muts.split('_')
        for i in mt:
            aa = i.strip()
            if (aa not in df_ratios.keys()):
                df_ratios[aa] = []
            df_ratios[aa].append(ratio)


    def setup_frame(self, df, substrate):
        df_values = {}
        for mut, ratio in zip(df[self.mutations], df[substrate]):
            self.set_values(mut, ratio, df_values)
        return df_values


    # generate csv file
    def get_csv_values(self, df, value, lst):
        for aa in df.keys():
            for val in df[aa]:
                tmp_ = aa + "," + str(val)
                lst.append(tmp_)
        return lst


    def setup_csv_file(self, substrate, df):
        tmpcsv = []
        value = substrate
        tmpcsv = self.get_csv_values(df, value, tmpcsv)
        with open(substrate + "_aa_act.csv", 'w') as f:
            f.write("AA,Ratio\n")
            for line in tmpcsv:
                f.write(line + "\n")
        rt_ = pd.read_csv(substrate + "_aa_act.csv")
        rt_ = rt_.dropna(axis=0, how='any')
        return rt_


    def bootstrap(self, df, aa, df_bootstrap, n=100):
        from random import randint
        for i in range(n):
            rd_ = randint(1, len(df[df['AA'] == aa]))
            if (aa not in df_bootstrap.keys()):
                df_bootstrap[aa] = []
            df_bootstrap[aa].append(df[df['AA'] == aa]["Ratio"].values[rd_ - 1])
        return df_bootstrap


    def get_bootstrap_frame(self, df):
        df_bootstrap = {}
        # bootstrap of values using N=10000
        for i in list(set(df["AA"])):
            df_bootstrap = self.bootstrap(df, i, df_bootstrap,self.nr_bootstraps)
        df_boot = pd.DataFrame.from_dict(df_bootstrap, orient='index').transpose()
        df_boot = df_boot.melt()
        df_boot.rename(columns={'variable': 'AA', 'value': 'Ratio'}, inplace=True)
        try:
            df_boot["pos"] = df_boot["AA"].str[1:-1]
            df_boot["pos"] = df_boot["pos"].astype(int)
        except:
            print("")
        df_boot.sort_values("pos", ascending=True, inplace=True)
        return df_boot

    def get_plotly_histogram_per_AA(self, df, value1, value2, name_of_plot1="test1", name_of_plot2="test2"):
        trace0 = go.Box(
            x=df["AA_"+value1],
            y=df["Ratio_"+value1],
            name=name_of_plot1,
            text="Variants: " + df["AA_"+value1],
        )

        trace1 = go.Box(
            x=df["AA_"+value2],
            y=df["Ratio_"+value2],
            name=name_of_plot2,
            text="Variants: " + df["AA_"+value2],
        )

        data = [trace0, trace1]

        layout = go.Layout(title="Bootstrap per AA", xaxis=dict(  # title="Amino acid substitution",
            titlefont=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f')),
                           yaxis=dict(
                               title="Ratio of activity compared with WT",
                               titlefont=dict(
                                   family='Courier New, monospace',
                                   size=18,
                                   color='#7f7f7f')),
                           boxmode='group')

        fig = go.Figure(data=data, layout=layout)
        return fig

    def get_aas(self, df, keyvalue='AA', value1='median', cutoff=1):
        activity = df[df[value1] > cutoff]
        return activity[keyvalue]


    def get_analysis_per_aa(self, df, values):
        '''

        :param values: list of values assigned to optimization e.g. ['r1', 'r2']
        :return:
        '''
        dfs = []
        for i in values:
            df_values_1 = self.setup_frame(df, i)
            rt_ = self.setup_csv_file(i, df_values_1)
            dfs.append(self.get_bootstrap_frame(rt_))

        df1 = dfs[0].reset_index(drop=True)
        df2 = dfs[1].reset_index(drop=True)

        dfnew = pd.merge(right=df1, left=df2,suffixes=("_"+values[0],"_"+values[1]), left_index=True, right_index=True)
        df_result1 = dfnew.groupby(by="AA_"+values[0], as_index=False)["Ratio_"+values[0]].agg([max, min, np.mean, gmean, np.std, np.median]).reset_index()
        df_result2 = dfnew.groupby(by="AA_"+values[1], as_index=False)["Ratio_"+values[1]].agg([max, min, np.mean, gmean, np.std, np.median]).reset_index()

        act1 = self.get_aas(df_result1,"AA_"+values[0],'median',self.cutoff)
        act2 = self.get_aas(df_result2,"AA_"+values[1],'median',self.cutoff)
        aas = set(list(act1) + list(act2))

        print(str(values[0])+" of samples: ", len(list(act1)))
        print(str(values[1])+" of samples: ", len(list(act2)))
        print("Number of unique AAs: ", len(aas))

        with open(self.today+"_designable_aas.txt", 'w') as f:
            for i in aas:
                f.write(i + "\n")
        fig = self.get_plotly_histogram_per_AA( dfnew, values[0], values[1], name_of_plot1=values[0], name_of_plot2=values[1])
        return fig, aas
