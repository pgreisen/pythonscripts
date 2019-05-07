from __future__ import print_function
import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import statsmodels.api as sm
import plotlyplots as pp


class EnrichmentAnalysis:

    df: pd.DataFrame

    def __init__(self):
        self._wt = None
        self.df = None
        self.pseudo_counts = 0.5
        # The name of the column of the naive library
        self.c0 = "c_0"
        # time_points
        self.timepoints = []
        # name or index of wildtype or reference sequence
        self.wt = ""
        self.wt_index = None
        self.wt_counts = {}
        self.seq_id_column = "seq_id"
        self.wt_scores = {}
        self.wt_se = {}
        self.file = ""
        self.column_name = "c_"


    def set_and_compute_z_and_p_value(self):
        for i in self.timepoints:
            self.df["Z_"+str(i) ] = (self.df["score_" + str(i)] - self.wt_scores["wt_score_"+str(i)]) / np.sqrt(self.df["SE_"+str(i)] ** 2 + self.wt_se["SE_"+str(i)] ** 2)
            self.df['pvalue_raw_'+str(i)] = 2 * scipy.stats.norm.sf(self.df["Z_"+str(i)])


    def set_wt_index(self):
        self.wt_index = self.df[self.df[self.seq_id_column] == self.wt].index


    def set_wt_counts(self):
        for i in self.timepoints:
            wt_ = "wt_" + str(i)
            self.wt_counts[wt_] = self.df.iloc[self.wt_index, self.df.columns.get_loc( "c_"+str(i))].values[0]

    def set_wt_logscores(self):
        for i in self.timepoints:
            wt_ = "wt_score_" + str(i)
            self.wt_scores[wt_] = self.df.iloc[self.wt_index, self.df.columns.get_loc("score_"+str(i))].values[0]

    def set_wt_se(self):
        for i in self.timepoints:
            wt_ = "SE_" + str(i)
            self.wt_se[wt_] = self.df.iloc[self.wt_index, self.df.columns.get_loc("SE_"+str(i))].values[0]


    def set_frequencies(self):
        '''
        :param df: dataframe with counts
        :param timepoints: list of column names for each time point
        :return: set the frequencies in the dataframe
        '''
        for i in self.timepoints:
            self.df["Freq_"+str(i)] = self.df[self.column_name+str(i)] / self.df[self.column_name+str(i)].sum()

    def compute_variant_ratio(self):
        '''

        :param self:
        :return:
        '''
        for i in self.timepoints:
            self.df["ratio_" + str(i)] = (self.df[self.column_name+str(i)] * self.wt_counts["wt_0"]) / (self.df[self.c0] * self.wt_counts["wt_" + str(i)])

    def compute_variant_log2_ratio(self):
        for i in self.timepoints:
            self.df["score_" + str(i)] = ((self.df[self.column_name+str(i)] + self.pseudo_counts) * (self.wt_counts["wt_0"] + self.pseudo_counts)) / \
                                         ((self.df[self.c0] + self.pseudo_counts) * (self.wt_counts["wt_" + str(i)] + self.pseudo_counts))

    def compute_variance_standard_errors(self):
        for i in self.timepoints:
            shared_variance =  1/ (self.wt_counts["wt_" + str(i)] + self.pseudo_counts) + 1/ (self.wt_counts["wt_0"] + self.pseudo_counts)
            self.df["Variance_" + str(i)] = np.sum(1. / (self.df[['c_0', 'c_'+str(i)]].values + self.pseudo_counts),axis=1) + shared_variance
            self.df["SE_" + str(i)] = np.sqrt(self.df["Variance_" + str(i)])



    def calc_regression(self, label):
        """
        Calculate least squares regression for *label*. If *weighted* is ``True``, calculates weighted least squares; else ordinary least squares.

        Regression results are stored in ``'/main/label/scores'``

        """
        if self.scoring_method == "WLS":
            for data in self.store.select_as_multiple(["/main/{}/log_ratios".format(label), "/main/{}/weights".format(label)], chunksize=self.chunksize):
                logging.info("Calculating weighted least squares for chunk {} ({} rows)".format(chunk, len(data.index)), extra={'oname' : self.name})
                result = data.apply(regression_apply, args=[self.timepoints, True], axis="columns")
                self.store.append("/main/{}/scores".format(label), result, min_itemsize={"index" : longest})
                chunk += 1
        elif self.scoring_method == "OLS":
            for data in self.store.select("/main/{}/log_ratios".format(label), chunksize=self.chunksize):
                logging.info("Calculating ordinary least squares for chunk {} ({} rows)".format(chunk, len(data.index)), extra={'oname' : self.name})
                result = data.apply(regression_apply, args=[self.timepoints, False], axis="columns")
                self.store.append("/main/{}/scores".format(label), result, min_itemsize={"index" : longest})
                chunk += 1
        else:
            raise ValueError('Invalid regression scoring method "{}" [{}]'.format(self.scoring_method, self.name))

        # need to read from the file, calculate percentiles, and rewrite it
        logging.info("Calculating slope standard error percentiles ({})".format(label), extra={'oname' : self.name})
        data = self.store['/main/{}/scores'.format(label)]
        data['score'] = data['slope']
        data['SE'] = data['SE_slope']
        data['SE_pctile'] = [stats.percentileofscore(data['SE'], x, "weak") for x in data['SE']]
        data = data[['score', 'SE', 'SE_pctile', 'slope', 'intercept', 'SE_slope', 't', 'pvalue_raw']] # reorder columns
        self.store.put("/main/{}/scores".format(label), data, format="table", data_columns=data.columns)


    def wt_plot(self, pdf):
        """
        Create a plot of the linear fit of the wild type variant.

        *pdf* is an open PdfPages instance.

        Only created for selections that use WLS or OLS scoring and have a wild type specified. 
        Uses :py:func:`~plots.fit_axes` for the plotting.
        """
        logging.info("Creating wild type fit plot", extra={'oname' : self.name})

        # get the data and calculate log ratios
        if "variants" in self.labels:
            wt_label = "variants"
        elif "identifiers" in self.labels:
            wt_label = "identifiers"
        data = self.store.select("/main/{}/counts".format(wt_label), where='index = "{}"'.format(WILD_TYPE_VARIANT)).ix[0]
        sums = self.store['/main/{}/counts'.format(wt_label)].sum(axis="index")  # sum of complete cases (N')
        yvalues = np.log(data + 0.5) - np.log(sums + 0.5)
        xvalues = [tp / float(max(self.timepoints)) for tp in self.timepoints]

        # fit the line
        X = sm.add_constant(xvalues) # fit intercept
        if self.scoring_method == "WLS":
            W =  1 / (1 / (data + 0.5) + 1 / (sums + 0.5))
            fit = sm.WLS(yvalues, X, weights=W).fit()
        elif self.scoring_method == "OLS":
            fit = sm.OLS(yvalues, X).fit()
        else:
            raise ValueError('Invalid regression scoring method "{}" [{}]'.format(self.scoring_method, self.name))
        intercept, slope = fit.params
        slope_se = fit.bse['x1']

        # make the plot
        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fit_axes(ax, xvalues, yvalues, slope, intercept, xlabels=self.timepoints)
        fit_axes_text(ax, cornertext="Slope {:3.2f}\nSE {:.1f}".format(slope, slope_se))
        ax.set_title("Wild Type Shape\n{}".format(self.name))
        ax.set_ylabel("Log Ratio (Complete Cases)")

        pdf.savefig(fig)
        plt.close(fig)


    def calc_outliers(self, label, minimum_components=4, log_chunksize=20000):
        """
        Test whether an element's individual components have significantly different 
        scores from the element. Results are stored in ``'/main/<label>/outliers'``.

        Args:
            label (str): label for the component (``'variants'`` or ``'barcodes'``)

            minimum_components (int): minimum number of componenents required for any statistics to be calculated

            log_chunksize (int): will output a log message every *n* rows

        """
        if self.check_store("/main/{}/outliers".format(label)):
            return

        if label == "variants":
            label2 = "synonymous"
        elif label == "barcodes":
            if self.is_barcodevariant():
                label2 = "variants"
            elif self.is_barcodeid():
                label2 = "identifiers"
            else:
                # this should never happen
                raise AttributeError("Failed to identify parent label when calculating barcode outliers [{}]".format(self.name))
        else:
            raise KeyError("Invalid label '{}' for calc_outliers [{}]".format(label,  self.name))

        logging.info("Identifying outliers ({}-{})".format(label, label2), extra={'oname' : self.name})
        
        
        logging.info("Mapping {} to {}".format(label, label2), extra={'oname' : self.name})
        if label == "variants":
            mapping = self.synonymous_variants()
        elif label == "barcodes":
            mapping = self.barcodemap_mapping()
        else:
            raise KeyError("Invalid label '{}' for calc_outliers [{}]".format(label,  self.name))

        # get the scores
        df1 = self.store.select("/main/{}/scores".format(label), "columns=['score', 'SE']")
        df2 = self.store.select("/main/{}/scores".format(label2), "columns=['score', 'SE']")

        # pre-calculate variances
        df1['var'] = df1['SE'] ** 2
        df1.drop('SE', axis=1, inplace=True)
        df2['var'] = df2['SE'] ** 2
        df2.drop('SE', axis=1, inplace=True)

        # set up empty results DF
        result_df = pd.DataFrame(np.nan, index=df1.index, columns=['z', 'pvalue_raw', 'parent'])

        # because this step can be slow, output chunk-style logging messages
        # pre-calculate the lengths for the log messages
        log_chunksize_list = [log_chunksize] * (len(df2) / log_chunksize) + [len(df2) % log_chunksize]
        log_chunk = 1

        for i, x in enumerate(df2.index):
            if i % log_chunksize == 0:
                logging.info("Calculating outlier p-values for chunk {} ({} rows) ({}-{})".format(log_chunk, log_chunksize_list[log_chunk - 1], label, label2), extra={'oname' : self.name})
                log_chunk += 1
            try:
                components = df1.loc[mapping[x]].dropna(axis="index", how="all")
            except KeyError:
                # none of the components were in the index
                continue
            if len(components.index) >= minimum_components:
                for c in components.index:
                    zvalue = np.absolute(df2.loc[x, 'score'] - df1.loc[c, 'score']) / np.sqrt(df2.loc[x, 'var'] + df1.loc[c, 'var'])
                    result_df.loc[c, 'z'] = zvalue
                    result_df.loc[c, 'pvalue_raw'] = 2 * stats.norm.sf(zvalue)
                    result_df.loc[c, 'parent'] = x
        if WILD_TYPE_VARIANT in result_df.index:
            result_df.loc[WILD_TYPE_VARIANT, 'z'] = np.nan
            result_df.loc[WILD_TYPE_VARIANT, 'pvalue_raw'] = np.nan
        result_df['z'] = result_df['z'].astype(float)
        result_df['pvalue_raw'] = result_df['pvalue_raw'].astype(float)

        self.store.put("/main/{}/outliers".format(label), result_df, format="table", data_columns=result_df.columns)

    def regression_apply(self, timepoints, name="score", weighted=False):

        """
        :py:meth:`pandas.DataFrame.apply` apply function for calculating
        enrichment using linear regression. If *weighted* is ``True`` perform
        weighted least squares; else perform ordinary least squares.

        Weights for weighted least squares are included in *row*.

        Returns a :py:class:`pandas.Series` containing regression coefficients,
        residuals, and statistics.
        """
        # retrieve log ratios from the row
        y = self.df[[name + '_{}'.format(t) for t in timepoints]]

        new_cols = ['OLS_intercept', 'OLS_slope', 'OLS_SE_slope', 'OLS_t', 'OLS_pvalue_raw'] + ['OLS_e_{}'.format(t) for t in timepoints]
        self.df = pd.concat([self.df, pd.DataFrame(columns=new_cols)],sort=False)
        # re-scale the x's to fall within [0, 1]
        xvalues = [x / float(max(timepoints)) for x in timepoints]

        # perform the fit
        X = sm.add_constant(xvalues)  # fit intercept
        if weighted:
            W = self.df[['W_{}'.format(t) for t in timepoints]]
            fit = sm.WLS(y, X, weights=W).fit()
        else:
            for index, row in y.iterrows():
                # 0: index
                # 1: values of df
                model = sm.OLS(np.array(row), X)
                results = model.fit()
                values = np.concatenate([results.params, [results.bse[0], results.tvalues[0], results.pvalues[0]], results.resid])

                for k,m in zip(new_cols, values):
                    self.df.iloc[index, self.df.columns.get_loc(k)] = m


    def main(self):
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        logger.info('Starting to populate dataframe')

        parser = argparse.ArgumentParser(description="Run NGS analysis.")
        parser.add_argument('-p', dest='pseudo_counts', help='Pseudo count add to low copy numbers', default=0.5, type=float)
        parser.add_argument('-t', dest='timepoints', help='Number of cycle for the NGS analysis using the format 1,2,3', nargs = '+', type = int )
        parser.add_argument('-f', dest='file',help='csv file with sequences and counts. Assumes counts are in columns like c_0 for naive library')
        parser.add_argument('-n', dest='wt', help="Identifier to identify wild type or reference sequence")
        parser.add_argument('-i', dest='seq_id_column', help="Name of column with IDs to identify wildtype or reference sequence", default="seq_id",type=str)


        args_dict = vars( parser.parse_args() )

        for item in args_dict:
            setattr(self, item, args_dict[item])

        # read input file
        # TODO : add assert statements to make sure data correctly read
        self.df = pd.read_csv(self.file)

        # set wild type id
        self.set_wt_index()

        # set counts for reference
        self.set_wt_counts()
        # compute frequencies for counts
        # this is not used at the moment
        self.set_frequencies()
        ##
        self.compute_variant_ratio()

        # Log2 of ratios with pseudo count
        self.compute_variant_log2_ratio()
        # Compute variance and standard error
        self.compute_variance_standard_errors()

        self.set_wt_logscores()
        self.set_wt_se()

        self.set_and_compute_z_and_p_value()


        logger.debug('WT index: %s', self.wt_index)
        logger.info('Set the counts of WT %s', self.wt_counts)
        # update records here
        logger.info('Finish populating dataframe')

        #print(self.df[['SE_1','SE_2','SE_3']])
        pp.barplot(self.df,['SE_1','SE_2','SE_3'],"Test")

        a = self.regression_apply( self.timepoints)
        self.df.to_excel("debug.xlsx", index=False)
        print(a)


if __name__ == "__main__":
        run = EnrichmentAnalysis()
        run.main()