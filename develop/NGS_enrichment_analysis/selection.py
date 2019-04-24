#  Copyright 2016-2017 Alan F Rubin
#
#  This file is part of Enrich2.
#
#  Enrich2 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Enrich2 is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Enrich2.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
import os
import re
import math
import time
import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as stats
import statsmodels.api as from .sfmap import sfmap_plot
from .plots import fit_axes, fit_axes_text, volcano_plot, configure_axes, plot_colors, weights_plot
from .constants import WILD_TYPE_VARIANT, SYNONYMOUS_VARIANT, AA_CODES
from .variant import protein_variant
from .dataframe import singleton_dataframe


class EnrichmentAnalysis:

    """
    Class for a single selection replicate, consisting of multiple 
    timepoints. This class coordinates :py:class:`~seqlib.seqlib.SeqLib` 
    objects.
    """
    

    def __init__(self):
        self._wt = None
        self.df = None
        self.pseudo_counts = 0.5
        # The name of the column of the naiive library
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


    def regression_apply(self, df, timepoints, name="score", weighted=False):
        """
        :py:meth:`pandas.DataFrame.apply` apply function for calculating
        enrichment using linear regression. If *weighted* is ``True`` perform
        weighted least squares; else perform ordinary least squares.

        Weights for weighted least squares are included in *row*.

        Returns a :py:class:`pandas.Series` containing regression coefficients,
        residuals, and statistics.
        """
        # retrieve log ratios from the row
        y = df[[name+'_{}'.format(t) for t in timepoints]]

        # re-scale the x's to fall within [0, 1]
        xvalues = [x / float(max(timepoints)) for x in timepoints]

        # perform the fit
        X = sm.add_constant(xvalues)  # fit intercept
        if weighted:
            W = df[['W_{}'.format(t) for t in timepoints]]
            fit = sm.WLS(y, X, weights=W).fit()
        else:
            fit = sm.OLS(y, X).fit()

        # re-format as a data frame row
        values = np.concatenate([fit.params, [fit.bse['x1'], fit.tvalues['x1'],
                                              fit.pvalues['x1']], fit.resid])
        index = ['intercept', 'slope', 'SE_slope', 't', 'pvalue_raw'] + \
                ['e_{}'.format(t) for t in timepoints]
        return pd.Series(data=values, index=index)


    def set_and_compute_z_and_p_value(self,base_score, base_se, col_val="log_ratio_3", col_se="se_ratio_3"):
        for i in self.timepoints:
            self.df["Z_"+i] = (self.df["score_" + str(i)] - self.wt_scores["wt_score_"+str(i)]) / np.sqrt(self.df["SE_"+str(i)] ** 2 + self.wt_se["SE_"+str(i))] ** 2)
            self.df['pvalue_raw_'+i] = 2 * scipy.stats.norm.sf(self.df["Z_"+str(i)])


    def set_wt_index(self):
        self.wt_index = self.df[self.df[self.seq_id_column] == self.wt].index


    def set_wt_counts(self):
        for i in self.timepoints:
            wt_ = "wt_" + str(i)
            self.wt_counts[wt_] = self.df.iloc[self.wt_index, i].value


    def set_wt_logscores(self):
        for i in self.timepoints:
            wt_ = "wt_score_" + str(i)
            self.wt_scores[wt_] = self.df.iloc[self.wt_index, self.df.columns.get_loc("score_"+str(i))].value


    def set_wt_se(self):
        for i in self.timepoints:
            wt_ = "wt_se_" + str(i)
            self.wt_se[wt_] = self.df.iloc[self.wt_index, self.df.columns.get_loc("SE_"+str(i))].value


    def set_frequencies(self,df,timepoints):
        '''
        :param df: dataframe with counts
        :param timepoints: list of column names for each time point
        :return: set the frequencies in the dataframe
        '''
        for i in self.timepoints:
            self.df["Freq_"+str(i)] = self.df[i] / self.df[i].sum()

    def compute_variant_ratio(self):
        '''

        :param self:
        :return:
        '''
        for i in self.timepoints:
            self.df["ratio_" + str(i)] = (self.df[i] * self.wt_counts["wt_c0"]) / (self.df[self.c0] * self.wt_counts["wt_" + str(i)])

    def compute_variant_log2_ratio(self):
        for i in self.timepoints:
            self.df["score_" + str(i)] = ((self.df[i] + self.pseudo_counts) * (self.wt_counts["wt_c0"] + self.pseudo_counts)) / \
                                         ((self.df[self.c0] + self.pseudo_counts) * (self.wt_counts["wt_" + str(i)] + self.pseudo_counts))

    def compute_variance_standard_errors(self):
        for i in self.timepoints:
            shared_variance =  1/ (self.wt_counts["wt_" + str(i)] + self.pseudo_counts) + 1/ (self.wt_counts["wt_c0"] + self.pseudo_counts)
            self.df["Variance_" + str(i)] = np.sum(1. / (df[['c_0', i]].values + self.pseudo_counts),axis=1) + shared_variance
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


    def se_pctile_plot(self, label, pdf):
        """
        Create plots of the linear fit of 21 selected variants, evenly distributed based on their standard error percentiles (0, 5, 10, ..., 100).

        *label* is the data label (barcode, variant, etc.)

        *pdf* is an open PdfPages instance.

        Uses :py:func:`~plots.fit_axes` for the plotting.
        """
        logging.info("Creating representative fit plots ({})".format(label), extra={'oname' : self.name})

        se_data = self.store.select("/main/{}/scores".format(label), where="columns in ['slope', 'intercept', 'SE_pctile'] & index!='{}' & index!='{}'".format(WILD_TYPE_VARIANT, SYNONYMOUS_VARIANT))
        se_data.sort_values("SE_pctile", inplace=True)
        se_data.dropna(axis="index", how="any", inplace=True)

        indices = np.linspace(0, len(se_data.index) - 1, 21).astype("int")
        se_data = se_data.ix[indices]

        # retrieves the whole DF because one case was hanging when trying to use select
        # totally unexplained, should fix later
        #ratio_data = self.store.select("/main/{}/log_ratios".format(label), "index=se_data.index")
        ratio_data = self.store.select("/main/{}/log_ratios".format(label)).loc[se_data.index]
        
        fig, axarr = plt.subplots(7, 3, sharex=True, sharey=True)
        fig.subplots_adjust(hspace=0, wspace=0) # eliminate white space between the subplots
        fig.set_size_inches((10, 17))
        fig.suptitle("Representative {} Fits\n{} ({})".format(self.scoring_method, self.name, label.title()))
        fig.subplots_adjust(top=0.958)           # eliminate white space after the title

        # turn off all tick labels
        plt.setp([ax.get_xticklabels() for ax in axarr.reshape(-1)], visible=False)
        plt.setp([ax.get_yticklabels() for ax in axarr.reshape(-1)], visible=False)

        # tick labels back on for some plots
        plt.setp([ax.get_xticklabels() for ax in (axarr[-1,0], axarr[-1,2])], visible=True)
        plt.setp([ax.get_yticklabels() for ax in (axarr[0,0], axarr[2,0], axarr[4,0], axarr[6,0])], visible=True)
        plt.setp([ax.get_yticklabels() for ax in (axarr[1,-1], axarr[3,-1], axarr[5,-1])], visible=True)
        plt.setp([ax.yaxis for ax in (axarr[1,-1], axarr[3,-1], axarr[5,-1])], ticks_position="right")

        # create the fit plots and add text to the individual plots
        for i, ax in enumerate(axarr.reshape(-1)):
            index = se_data.index[i]
            fit_axes(ax, xvalues=[x / float(max(self.timepoints)) for x in self.timepoints], yvalues=ratio_data.loc[index], slope=se_data.loc[index, "slope"], intercept=se_data.loc[index, "intercept"], xlabels=self.timepoints)
            fit_axes_text(ax, cornertext="Slope {:3.2f}\nSE Pctile {:.1f}".format(se_data.loc[index, "slope"], se_data.loc[index, "SE_pctile"]), centertext=index)

        # turn off the subplot axis labels
        [ax.set_xlabel("") for ax in axarr.reshape(-1)]
        [ax.set_ylabel("") for ax in axarr.reshape(-1)]

        # add x and y labels in the middle
        axarr[-1, 1].set_xlabel("Time Point")
        axarr[3, 0].set_ylabel("Log Ratio")

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


    def timepoint_counts_plot(self, label, pdf):
        """
        Create barplot of the number of items counted for each time point.

        *label* is the data label (barcode, variant, etc.)

        *pdf* is an open PdfPages instance.
        """
        logging.info("Creating time point count plots ({})".format(label), extra={'oname' : self.name})

        counts = self.store['/main/{}/counts'.format(label)].sum(axis="index")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        configure_axes(ax)

        xpos = np.arange(len(self.timepoints))
        width = 0.8
        ax.bar(xpos, counts, width, color=plot_colors['bright1'])
        ax.set_title("Total {}\n{}".format(label.title(), self.name))
        ax.set_ylabel("Count")
        ax.set_xlabel("Timepoint")
        ax.set_xticks(xpos + width / 2.)
        ax.set_xticklabels(self.timepoints)

        pdf.savefig(fig)
        plt.close(fig)


    def volcano_plot(self, label, pdf, colors="YlGnBu_r", log_bins=True):
        """
        Create a volcano plot (p-value vs. functional score).

        *label* is the data label (barcode, variant, etc.)

        *pdf* is an open PdfPages instance.

        The p-values used are the regression p-values (p-value of non-zero slope). Due to the large number of points, we use a hexbin plot showing the density instead of a scatter plot.
        """
        logging.info("Creating volcano plot ({})".format(label), extra={'oname' : self.name})

        # get the data
        data = self.store.select("/main/{}/scores".format(label), "columns=['score', 'pvalue_raw']")
        volcano_plot(data, pdf, title="{} ({})".format(self.name, label.title()), colors=colors, log_bins=log_bins)


    def make_plots(self):
        """
        Create plots for this entity.

        This function handles opening and closing the various PDF files for multi-page plots, as well as plotting similar plots for different data labels.
        """
        if self.plots_requested:
            logging.info("Creating plots", extra={'oname' : self.name})

            # counts per time point
            pdf = PdfPages(os.path.join(self.plot_dir, "timepoint_counts.pdf"))
            for label in self.labels:
                self.timepoint_counts_plot(label, pdf)
            pdf.close()

            # wild type shape
            if self.logr_method == "wt" and self.scoring_method in ("WLS", "OLS"):
                pdf = PdfPages(os.path.join(self.plot_dir, "wt_shape.pdf"))
                self.wt_plot(pdf)
                pdf.close()

            # regression weights
            if self.scoring_method == "WLS":
                pdf = PdfPages(os.path.join(self.plot_dir, "regression_weights.pdf"))
                for label in self.labels:
                    weights_plot(self, label, pdf)
                pdf.close()

            # linear fits by standard error percentile
            if self.scoring_method in ("WLS", "OLS"):
                pdf = PdfPages(os.path.join(self.plot_dir, "se_pctile.pdf"))
                for label in self.labels:
                    self.se_pctile_plot(label, pdf)
                pdf.close()

            # volcano plots
            #if self.scoring_method in ("WLS", "OLS", "ratios"):
            if self.scoring_method in ("WLS", "OLS") and "variants" in self.labels:
                pdf = PdfPages(os.path.join(self.plot_dir, "volcano.pdf"))
                for label in self.labels:
                    self.volcano_plot(label, pdf, log_bins=True)
                pdf.close()

            # library diversity for each time point (amino acid)
            if "synonymous" in self.labels:
                pdf = PdfPages(os.path.join(self.plot_dir, "diversity_aa.pdf"))
                for tp in self.timepoints:
                    self.sfmap_wrapper(cname="c_{}".format(tp), pdf=pdf, coding=True, log10=True)
                pdf.close()

            # library diversity for each time point (nucleotide)
            if "variants" in self.labels:
                pdf = PdfPages(os.path.join(self.plot_dir, "diversity_nt.pdf"))
                for tp in self.timepoints:
                    self.sfmap_wrapper(cname="c_{}".format(tp), pdf=pdf, coding=False, log10=True)
                pdf.close()

            # sequence-function maps
            if self.scoring_method != "counts":
                if "synonymous" in self.labels:
                    pdf = PdfPages(os.path.join(self.plot_dir, "sequence_function_map_aa.pdf"))
                    self.sfmap_wrapper(cname="score", pdf=pdf, coding=True)
                    pdf.close()
                if "variants" in self.labels:
                    pdf = PdfPages(os.path.join(self.plot_dir, "sequence_function_map_nt.pdf"))
                    self.sfmap_wrapper(cname="score", pdf=pdf, coding=False)
                    pdf.close()

        # SeqLib plots
        for lib in self.children:
            lib.make_plots()


    def write_tsv(self):
        """
        Write each table from the store to its own tab-separated file.

        Files are written to a ``tsv`` directory in the default output location. 
        File names are the HDF5 key with ``'_'`` substituted for ``'/'``.
        """
        if self.tsv_requested:
            logging.info("Generating tab-separated output files", extra={'oname' : self.name})
            for k in self.store.keys():
                self.write_table_tsv(k)
        for lib in self.children:
            lib.write_tsv()


    def synonymous_variants(self):
        """
        Populate and return a dictionary mapping synonymous variants to the 
        list of associated variants in ``/main/variants/counts``.
        """
        mapping = dict()
        try:
            variants = self.store.select_column("/main/variants/counts", "index")
        except KeyError:
            raise KeyError("No variant counts found [{}]".format(self.name))
        for v in variants:
            pv = protein_variant(v)
            try:
                mapping[pv].append(v)
            except KeyError:
                mapping[pv] = [v]
        return mapping


    def sfmap_wrapper(self, cname, pdf, coding, log10=False):
        """
        Create a sequence function map for either scores or library diversity.

        Uses :py:func:`~sfmap.sfmap_plot` for the plotting.
        """
        plot_options = self.get_root().plot_options
        
        if cname.startswith("c_"):
            counts = True
        elif cname == "score":
            counts = False
        else:
            raise ValueError("Invalid sequence-function map data column.")

        if coding:
            label = "amino acid"
        else:
            label = "nucleotide"

        if counts:
            logging.info("Creating diversity map ({})".format(label), extra={'oname' : self.name})
        else:
            logging.info("Creating sequence-function map ({})".format(label), extra={'oname' : self.name})

        # build the data frame name and get the data
        df_name = "/main/"
        if coding:
            df_name += "synonymous/"
        else:
            df_name += "variants/"
        if counts:
            df_name += "counts_unfiltered"
        else:
            df_name += "scores"
        if plot_options is not None:
            data, wtseq = singleton_dataframe(self.store[df_name][cname], self.wt, coding=coding, aa_list=plot_options['aa_list'])
        else:
            data, wtseq = singleton_dataframe(self.store[df_name][cname], self.wt, coding=coding)
        if counts:
            data_se = None
        else:
            if plot_options is not None:
                data_se, _ = singleton_dataframe(self.store[df_name]["SE"], self.wt, coding=coding, aa_list=plot_options['aa_list'])
            else:
                data_se, _ = singleton_dataframe(self.store[df_name]["SE"], self.wt, coding=coding)


        # format the title
        if coding:
            title = "Amino Acid"
        else:
            title = "Nucleotide"
        if counts:
            title += " Diversity Map\n{} (Time {})".format(self.name, cname[2:]) # trim the "c_"
        else:
            if self.scoring_method in ("WLS", "OLS"):
                title += " Sequence-Function Map\n{} ({} Slope)".format(self.name, self.scoring_method)
            elif self.scoring_method == "ratios":
                title += " Sequence-Function Map\n{} ({})".format(self.name, "Enrich2 Ratio")
            elif self.scoring_method == "simple":
                title += " Sequence-Function Map\n{} ({})".format(self.name, "Simplified Ratio")
            else:
                raise ValueError("Invalid scoring method", self.name)

        if counts and log10:
            style = "logcounts"
        elif counts:
            style = "counts"
        else:
            style = "scores"

        if plot_options is not None:
            sfmap_plot(df=data, pdf=pdf, style=style, df_se=data_se,
                       dimensions="tall", wt=wtseq, title=title,
                       aa_list=plot_options['aa_list'],
                       aa_label_groups=plot_options['aa_label_groups'])
        else:
            sfmap_plot(df=data, pdf=pdf, style=style, df_se=data_se,
                       dimensions="tall", wt=wtseq, title=title)



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




    def volcano_plot(self, label, pdf, colors="YlGnBu_r", log_bins=True):
        """
        Create a volcano plot (p-value vs. functional score).

        *label* is the data label (barcode, variant, etc.)

        *pdf* is an open PdfPages instance.

        The p-values used are the regression p-values (p-value of non-zero slope). Due to the large number of points,
        we use a hexbin plot showing the density instead of a scatter plot.

        """
        logging.info("Creating volcano plot ({})".format(label), extra={'oname' : self.name})

        # get the data
        data = self.store.select("/main/{}/scores".format(label), "columns=['score', 'pvalue_raw']")
        volcano_plot(data, pdf, title="{} ({})".format(self.name, label.title()), colors=colors, log_bins=log_bins)

    def regression_apply(self, row, timepoints, weighted):
        """
        :py:meth:`pandas.DataFrame.apply` apply function for calculating
        enrichment using linear regression. If *weighted* is ``True`` perform
        weighted least squares; else perform ordinary least squares.

        Weights for weighted least squares are included in *row*.

        Returns a :py:class:`pandas.Series` containing regression coefficients,
        residuals, and statistics.
        """
        # retrieve log ratios from the row
        y = row[['L_{}'.format(t) for t in timepoints]]

        # re-scale the x's to fall within [0, 1]
        xvalues = [x / float(max(timepoints)) for x in timepoints]

        # perform the fit
        X = sm.add_constant(xvalues)  # fit intercept
        if weighted:
            W = row[['W_{}'.format(t) for t in timepoints]]
            fit = sm.WLS(y, X, weights=W).fit()
        else:
            fit = sm.OLS(y, X).fit()

        # re-format as a data frame row
        values = np.concatenate([fit.params, [fit.bse['x1'], fit.tvalues['x1'],
                                              fit.pvalues['x1']], fit.resid])
        index = ['intercept', 'slope', 'SE_slope', 't', 'pvalue_raw'] + \
                ['e_{}'.format(t) for t in timepoints]
        return pd.Series(data=values, index=index)

    def calc_weights(self, label):
        """
        Calculate the regression weights (1 / variance).
        """

        # perform operations on the numpy values of the data frame for easier broadcasting
        #variances = 1.0 / (variances[c_n].values + 0.5) + 1.0 / (variances[['c_0']].values + 0.5)
        variances = 1.0 / (variances[c_n].values + 0.5)
        if self.logr_method == "wt":
            if "variants" in self.labels:
                wt_label = "variants"
            elif "identifiers" in self.labels:
                wt_label = "identifiers"
            else:
                raise ValueError('Failed to use wild type log ratio method, suitable data table not present [{}]'.format(self.name))
            wt_counts = self.store.select("/main/{}/counts".format(wt_label), "columns=c_n & index='{}'".format(WILD_TYPE_VARIANT))
            if len(wt_counts) == 0: # wild type not found
                raise ValueError('Failed to use wild type log ratio method, wild type sequence not present [{}]'.format(self.name))
            variances = variances + 1.0 / (wt_counts.values + 0.5)
        elif self.logr_method == "complete":
            variances = variances + 1.0 / (self.store.select("/main/{}/counts".format(label), "columns=c_n").sum(axis="index").values + 0.5)
        elif self.logr_method == "full":
            variances = variances + 1.0 / (self.store.select("/main/{}/counts_unfiltered".format(label), "columns=c_n").sum(axis="index", skipna=True).values + 0.5)
        else:
            raise ValueError('Invalid log ratio method "{}" [{}]'.format(self.logr_method, self.name))
        variances = 1.0 / variances # weights are reciprocal of variances
        # make it a data frame again
        variances = pd.DataFrame(data=variances, index=index, columns=['W_{}'.format(x) for x in self.timepoints])
        self.store.put("/main/{}/weights".format(label), variances, format="table", data_columns=variances.columns)
