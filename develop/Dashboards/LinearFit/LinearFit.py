import pandas as pd
import seaborn as sns
import numpy as np
from scipy.stats import gmean
import matplotlib.pyplot as plt
import math
import datetime

from datetime import datetime as dt
import numpy as np
import pandas as pd
from ipywidgets import widgets
from plotly.subplots import make_subplots
from plotly.graph_objs import Scatter, Layout
import plotly
import plotly.offline as offline
import plotly.graph_objs as go
import plotly.plotly as py
import plotly.figure_factory as ff
plotly.offline.init_notebook_mode(connected=False)


class LinearFit:

    def __init__(self):
        self.date = str(datetime.date.today())
        # Set a cutoff to remove bad fits.
        # Alpha > 0;
        # R2 > 0.85;
        # Sort by R2, alpha, window of fit
        self.r2_cutoff = 0.85
        self.alpha_cutoff = 0
        ######## Parameters for the fit ########
        # slicing fits into smaller pieces
        self.minimum_window = [10,12,15,20]
        # number of top fits to export
        self.number_of_frame_to_export = 2
        # Threshold if bigger than this raise a flag for the fit
        self.threshold = 25.0
        self.inp_file = "2020-01-02_20191219_APassay_protease_challenge_AP-HT_variantsrepeat_project_KineticMeasurements.csv"
        # names for output file
        self.exp="elanco_20191219_APassay_protease_challenge_AP-HT_variantsrepeat"
        self.read = 2
        self.run="run_"+str(self.read)+"_"
        self.pth = "./"
        self.cols_to_keep = ['temperature', 'well', 'fluorescence', 'time_s']
        self.cols_to_summary = ['Exp','alpha','R2','CV','End','Range','Rank','Residual','Start','b','Attention']
        self.time_for_fit = 'time_min'
        self.dfs = None
        self.fitdata = None
        self.outputname = ""

    def get_dataframe(self):
        df = pd.read_csv(self.pth+self.inp_file)
        return df


    def get_filtered_read(self,df ):
        df = df[df['read'] == self.read]
        df = df[self.cols_to_keep]
        return df

    def get_reformatted_dataframe(self,df):
        # replace overflow with 0
        df['fluorescence']  = df['fluorescence'].replace('OVRFLW',0)
        df['fluorescence'] = df['fluorescence'].astype(float)
        return df

    def get_pivoted_table(self,df):
        df = df.pivot_table(index=['time_s', 'temperature'],columns=['well'])
        return df


    def get_order_of_output_columns(self,df_columns):
        order_of_output_columns = sorted(list(df_columns), key=lambda x: (x[0], x[1:-1]))
        return order_of_output_columns

    def get_minutes_column(self,df):
        # this has been changed to seconds and is converted into minutes here
        df["time_min"] = df['time_s'] / 60
        return df


    def get_linear_fit(self,time_spoint,y):
        # a : fitparameters ( slope/intersection)
        # residuals : residuals
        # rank : rank of matrix
        # sv : singular values
        # ct : conditioning threshold
        a, residuals,rank, sv,ct  = np.polyfit(time_spoint, y, 1, full=True)
        yfit = np.polyval(a,time_spoint)
        r2 = np.corrcoef(y,yfit)[0,1]**2
        return a, residuals, rank, sv, ct, time_spoint, y, yfit, r2



    def fit_exp(self, df, col, minimum_window,timevalue='time_min'):
        '''Returns dictionary with fitted values'''
        # Make sure that index starts from 0!
        number_of_obs = len(df.index)
        fits = {}
        for i in df.index:
            idx_length = i+minimum_window
            if(idx_length < number_of_obs):
                t = df[timevalue][i:idx_length]
                obs = df[col][i:idx_length].values
                key = "Interval_"+str(i)+"_"+str(idx_length)
                vals = self.get_linear_fit(t,obs)
                fits[key] = vals
        return fits


    def generate_dataframe(self,fits, col,range_of_fit, df):
        keys = list(fits.keys())
        # generate dataframe object
        cols = ['Exp','Start','End','alpha','b','Residual','R2','Rank','CV','Range']
        df_fit = pd.DataFrame(index=range(0,len(keys)),columns=cols)
        for i in range(0,len(keys)):
            tmp_interval = keys[i].split("_")
            start,end = tmp_interval[1],tmp_interval[2]
            tmp_tuple = fits[keys[i]]
            # Fit parameters
            df_fit.iloc[i,df_fit.columns.get_loc("Exp")] = col
            df_fit.iloc[i,df_fit.columns.get_loc("Start")] = df['time_s'][int(start)]
            df_fit.iloc[i,df_fit.columns.get_loc("End")] = df['time_s'][int(end)]
            df_fit.iloc[i,df_fit.columns.get_loc("alpha")] = tmp_tuple[0][0]
            df_fit.iloc[i,df_fit.columns.get_loc("b")] = tmp_tuple[0][1]
            df_fit.iloc[i,df_fit.columns.get_loc("R2")] = tmp_tuple[-1]
            df_fit.iloc[i,df_fit.columns.get_loc("Residual")] = tmp_tuple[1][0]
            df_fit.iloc[i,df_fit.columns.get_loc("CV")] = tmp_tuple[4]
            df_fit.iloc[i,df_fit.columns.get_loc("Rank")] = tmp_tuple[2]
            df_fit.iloc[i,df_fit.columns.get_loc("Range")] = range_of_fit
        return df_fit


    def get_formatted_dataframes(self, df):
        # list_of_dictionaries = []
        dataframes = {}
        length_of_dataframe = len(df.index)
        fitdata = {}
        # get fits for each column
        for i in df.columns[1:]:
            for j in self.minimum_window:
                fit_ = self.fit_exp(df,i,j)
                # need to fix this
                fitdata[i] = fit_
                tmp_key = str(i)+"_"+str(j)
                dataframes[tmp_key] =  self.generate_dataframe(fit_, i, j, df)

        # concatenate dataframe per experiment
        dfs = {}
        for i in df.columns[1:]:
            tmp = []
            for j in self.minimum_window:
                tmp.append(dataframes[i+"_"+str(j)])
            tmpdf = pd.concat(tmp)
            # Filter based on criteria
            tmpdf = tmpdf[(tmpdf["R2"] > self.r2_cutoff) & (tmpdf["alpha"] > self.alpha_cutoff)]
            dfs[i] = tmpdf.reset_index()

        # rename it to match earlier naming
        # dataframes = dfs
        return dfs, fitdata

    def get_fitted_data(self):
        return self.dfs, self.fitdata

    ##################
    ## Compute fits ##
    ##################

    def generate_empty_dataframe(self, col, dummy_string):
        # generate dataframe object
        cols = ['Exp', 'Start', 'End', 'alpha', 'b', 'Residual', 'R2', 'Rank', 'CV', 'Range', 'Attention']
        df_fit = pd.DataFrame(index=range(0, 1), columns=cols)
        i = 0
        # Fit parameters
        df_fit.iloc[i, df_fit.columns.get_loc("Exp")] = col
        df_fit.iloc[i, df_fit.columns.get_loc("Start")] = dummy_string
        df_fit.iloc[i, df_fit.columns.get_loc("End")] = dummy_string
        df_fit.iloc[i, df_fit.columns.get_loc("alpha")] = dummy_string
        df_fit.iloc[i, df_fit.columns.get_loc("b")] = dummy_string
        df_fit.iloc[i, df_fit.columns.get_loc("R2")] = dummy_string
        df_fit.iloc[i, df_fit.columns.get_loc("Residual")] = dummy_string
        df_fit.iloc[i, df_fit.columns.get_loc("CV")] = dummy_string
        df_fit.iloc[i, df_fit.columns.get_loc("Rank")] = dummy_string
        df_fit.iloc[i, df_fit.columns.get_loc("Range")] = dummy_string
        df_fit.iloc[i, df_fit.columns.get_loc("Attention")] = True
        return df_fit

    def write_fits_to_file(self, date, exp, run, cols_to_summary, minimum_window, dataframes,
                               number_of_frame_to_export,order_of_output_columns):
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        self.outputname = date + '_' + exp + '_' + run + '.xlsx'
        writer = pd.ExcelWriter( self.outputname, engine='xlsxwriter')
        # Write each dataframe to a different worksheet. you could write different string like above if you want
        dataframes['summary'][cols_to_summary].to_excel(writer, sheet_name='Summary', index=False)
        for i in order_of_output_columns[:-1]:
            if (i == 'summary'):
                continue
            if (i == 'time_s'):
                continue
            df_tofile = {}
            for j in minimum_window:
                dftmp = dataframes[i][dataframes[i]["Range"] == j]
                key = i + "_" + str(j)
                df_tofile[key] = dftmp.sort_values(by="R2", ascending=False).head(number_of_frame_to_export)
            dftmp = pd.concat(df_tofile)
            dftmp.sort_values(by="R2", ascending=False).to_excel(writer, sheet_name=i, index=False)
        # Close the Pandas Excel writer and output the Excel file.
        writer.save()

    def get_fit(self, order_of_output_columns, dataframes, minimum_window, number_of_frame_to_export):
        # order of columns
        dummy_string = "No fit possible"
        # get top slope
        df_tofile = []
        key = "summary"
        for i in order_of_output_columns:
            flag_ = False
            if (i == 'time_sPoints'):
                continue
            elif (i == 'time_s'):
                continue
            elif (i in dataframes.keys()):
                if (len(dataframes[i]) == 0):
                    df_tofile.append(self.generate_empty_dataframe(i, dummy_string))
                else:
                    # evaluate statistics here
                    tmp_ = []
                    for j in minimum_window:
                        dftmp = dataframes[i][dataframes[i]["Range"] == j]
                        tmp_.append(dataframes[i].sort_values(by="R2", ascending=False).head(number_of_frame_to_export))
                    flag_df = pd.concat(tmp_)
                    var_ = (flag_df['alpha'].max() - flag_df['alpha'].min()) / flag_df['alpha'].max() * 100
                    if (var_ > self.threshold):
                        flag_ = True
                    tmp_ = dataframes[i].sort_values(by="R2", ascending=False).head(1)
                    tmp_["Attention"] = flag_
                    df_tofile.append(tmp_)
            else:
                df_tofile.append(self.generate_empty_dataframe(i, dummy_string))

        dftmp = pd.concat(df_tofile)
        dataframes[key] = dftmp
        return dataframes

    def main(self):
        df = self.get_dataframe()
        # get read 1 e.g.
        df = self.get_filtered_read(df)
        # set overflow values to 0
        df = self.get_reformatted_dataframe(df)
        df = self.get_pivoted_table(df)

        # Remove the time column
        df.columns = df.columns.droplevel(0)
        df = df.rename_axis(None, axis=1)
        df = df.reset_index()
        df = df.drop(labels=['temperature'],axis=1)

        order_of_columns = self.get_order_of_output_columns(df.columns)
        df = self.get_minutes_column(df)
        df.fillna(value=0.0, inplace=True)
        self.dfs, self.fitdata = self.get_formatted_dataframes(df)

        # write results to excel file
        dataframe_w_fit = self.get_fit( order_of_columns, self.dfs, self.minimum_window, self.number_of_frame_to_export)
        self.write_fits_to_file(self.date, self.exp, self.run, self.cols_to_summary, self.minimum_window, dataframe_w_fit, self.number_of_frame_to_export,order_of_columns)



if __name__ == "__main__":
   run = LinearFit()
   run.main()