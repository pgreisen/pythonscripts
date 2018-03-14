import plotly
import plotly.offline as offline
import plotly.graph_objs as go
import plotly.plotly as py
import plotly.figure_factory as ff


class JupyterPlotly:


    def __init__(self):
        self.nv_min = 3
        self.nv_max = 6

    def get_scatter_plot_w_error_bars_index(self,dataframe, yvalue, std_err, spl_var, nnc, title, xlabel, ylabel,mut_label,*args):
        '''
        :param dataframe:
        :param xvalue:
        :param yvalue:
        :param std_err:
        :param spl_var:
        :param nnc:
        :param title:
        :param xlabel:
        :param ylabel:
        :param args:
        :return: plotly figure object
        '''
        data = []
        import types
        assert isinstance(args, types.TupleType)
        for i in args[0]:
            tmp_ = dataframe[(dataframe[spl_var] == i)]
            data.append(go.Scatter(y=tmp_[yvalue], \
 \
                                   error_y=dict(type='data', array=tmp_[std_err], visible=True), \
                                   mode='markers', text=tmp_[nnc] + "\n" + \
                                                        tmp_[mut_label],
                                   name="conc: " + str(i) + " nM", marker=dict(size=10))
                        )
        layout = go.Layout(title=title,
                           xaxis=dict(title=xlabel,
                                      titlefont=dict(
                                          family='Courier New, monospace',
                                          size=18,
                                          color='#7f7f7f'
                                      )
                                      ),
                           yaxis=dict(
                               title=ylabel,
                               titlefont=dict(
                                   family='Courier New, monospace',
                                   size=18,
                                   color='#7f7f7f')))

        return go.Figure(data=data, layout=layout)

    def get_scatter_plot_w_error_bars(self, dataframe, xvalue, yvalue, std_err, spl_var, nnc, title, xlabel, ylabel,
                                      mut_label, *args):
        '''
        :param dataframe:
        :param xvalue:
        :param yvalue:
        :param std_err:
        :param spl_var:
        :param nnc:
        :param title:
        :param xlabel:
        :param ylabel:
        :param args:
        :return: plotly figure object
        '''
        # CONC_ = u'Assay Conc (nM)'
        # xvalue = u'BiAb AC-SINS (%)'
        # yvalue = "gmean"
        # std_err = "std"
        # nnc = "NNC"
        # mut_label = "FX Mut (HC|LC)"
        # title = "AC-SINS versus TGT normalized"
        # xlabel = 'AC-SINS score(%)'
        # ylabel = 'Ratio (WT)'
        data = []
        import types
        assert isinstance(args, types.TupleType)
        for i in args[0]:
            tmp_ = dataframe[(dataframe[spl_var] == i)]
            data.append(go.Scatter(x=tmp_[xvalue],
                                   y=tmp_[yvalue], \
                                   error_y=dict(type='data', array=tmp_[std_err], visible=True), \
                                   mode='markers', text=tmp_[nnc] + "\n" + \
                                                        tmp_[mut_label],
                                   name="conc: " + str(i) + " nM", marker=dict(size=10))
                        )
        layout = go.Layout(title=title,
                           xaxis=dict(title=xlabel,
                                      titlefont=dict(
                                          family='Courier New, monospace',
                                          size=18,
                                          color='#7f7f7f'
                                      )
                                      ),
                           yaxis=dict(
                               title=ylabel,
                               titlefont=dict(
                                   family='Courier New, monospace',
                                   size=18,
                                   color='#7f7f7f')))

        return go.Figure(data=data, layout=layout)



    def get_scatter_plotly_size(self,dataframe,xlabel, ylabel,zlabel,ztextlabel,nnc, lc_label, hc_label,plottitle=""):
        '''
        :param dataframe:
        :param xlabel:
        :param ylabel:
        :param nnc:
        :param lc_label:
        :param hc_label:
        :param plottitle:
        :return:
        '''
        data = []
        # this formula is used to calculate the size of the bubble
        # sizeref = 2. * max(array of size values) / ( desired maximum marker size ** 2)
        size_ = self.get_scaled_values(list(dataframe[zlabel]))

        sizeref = 2 * max(size_) / (self.nv_max ** 2)

        data.append(go.Scatter(x=dataframe[xlabel], \
                               y=dataframe[ylabel], \
                               mode='markers',
                               text=dataframe[nnc] + "\n" + "LC: " \
                                    + dataframe[lc_label] + "\nHC: " + \
                                    dataframe[hc_label] + '\n'+ztextlabel+': ' + dataframe[zlabel].apply(str),
                               marker=dict(
                                   color=dataframe[zlabel],
                                   size=size_,
                                   sizeref=sizeref,
                                   showscale=True,
                                   colorbar=dict(title=ztextlabel)

                               )))
        layout = go.Layout(title=plottitle, xaxis=dict(title=xlabel,
                                                       titlefont=dict(
                                                           family='Courier New, monospace',
                                                           size=18,
                                                           color='#7f7f7f')),
                           yaxis=dict(
                               title=ylabel,
                               titlefont=dict(
                                   family='Courier New, monospace',
                                   size=18,
                                   color='#7f7f7f')))
        fig = go.Figure(data=data, layout=layout)
        return fig


    def get_scatter_plotly(self,dataframe,xlabel, ylabel, nnc, lc_label, hc_label,plottitle=""):
        '''

        :param dataframe:
        :param xlabel:
        :param ylabel:
        :param nnc:
        :param lc_label:
        :param hc_label:
        :param plottitle:
        :return:
        '''
        data = []
        data.append(go.Scatter(x=dataframe[xlabel],\
                               y=dataframe[ylabel],\
                               mode='markers',
                               text=dataframe[nnc] + "\n" + \
                                                    "LC: "\
                                    +dataframe[lc_label]+"\nHC: " + dataframe[hc_label],
                               marker=dict(size=10)))

        layout = go.Layout( title=plottitle,
                           xaxis=dict(title=xlabel,
                                      titlefont=dict(
                                          family='Courier New, monospace',
                                          size=18,
                                          color='#7f7f7f')),
                           yaxis=dict(
                               title=ylabel,
                               titlefont=dict(
                                   family='Courier New, monospace',
                                   size=18,
                                   color='#7f7f7f')))
        fig = go.Figure(data=data, layout=layout)
        return fig


    def get_scaled_values(self,listofvalues):
        ov_min = min(listofvalues)
        ov_max = max(listofvalues)
        tmp_ = []
        scale = (self.nv_max - self.nv_min) / (ov_max - ov_min)
        for i in listofvalues:
            nw_ = (i - ov_min) * scale + self.nv_min
            tmp_.append(nw_)
        return tmp_

    def get_scatter_plotly_muts(self, dataframe, xvalue, yvalue,spl_var, nnc, title, xlabel, ylabel,
                                      mut_label):
        '''
        :param dataframe:
        :param xvalue:
        :param yvalue:
        :param std_err:
        :param spl_var:
        :param nnc:
        :param title:
        :param xlabel:
        :param ylabel:
        :param args:
        :return: plotly figure object
        '''
        # CONC_ = u'Assay Conc (nM)'
        # xvalue = u'BiAb AC-SINS (%)'
        # yvalue = "gmean"
        # std_err = "std"
        # nnc = "NNC"
        # mut_label = "FX Mut (HC|LC)"
        # title = "AC-SINS versus TGT normalized"
        # xlabel = 'AC-SINS score(%)'
        # ylabel = 'Ratio (WT)'
        data = []

        data.append(go.Scatter(x=dataframe[xvalue],
                                   y=dataframe[yvalue], \
                                   mode='markers', text=dataframe[nnc] + "\n" + \
                                                        dataframe[mut_label],marker=dict(size=10))
                        )
        layout = go.Layout(title=title,
                           xaxis=dict(title=xlabel,
                                      titlefont=dict(
                                          family='Courier New, monospace',
                                          size=18,
                                          color='#7f7f7f'
                                      )
                                      ),
                           yaxis=dict(
                               title=ylabel,
                               titlefont=dict(
                                   family='Courier New, monospace',
                                   size=18,
                                   color='#7f7f7f')))

        return go.Figure(data=data, layout=layout)



    def get_scatter_plotly_z_heatmap(self, dataframe, xvalue, yvalue, zvalue, zlabel_text, spl_var, nnc,\
                                     lc_label, hc_label,title, xlabel, ylabel, *args):
        data = []
        import types
        assert isinstance(args, types.TupleType)
        for i in args[0]:
            tmp_ = dataframe[ (dataframe[spl_var] == i)]
            size_ = self.get_scaled_values(list(dataframe[zlabel]))
            sizeref = 2 * max(size_) / (self.nv_max ** 2)

            data.append(go.Scatter(x=tmp_[xlabel],y=tmp_[ylabel], \
                           mode='markers',
                           text=tmp_[nnc] + "\n" + "LC: " \
                                + tmp_[lc_label] + "\nHC: " + tmp_[hc_label] \
                                + '\n'+str(zlabel_text)+': ' + tmp_[zlabel].apply(str),
                           marker=dict(
                               color=tmp_[zlabel],
                               size=size_,
                               sizeref=sizeref,
                               showscale=True,
                               colorbar=dict(title=ztextlabel)

                           )))
            layout = go.Layout(title=plottitle, xaxis=dict(title=xlabel,
                                                   titlefont=dict(
                                                       family='Courier New, monospace',
                                                       size=18,
                                                       color='#7f7f7f')),
                       yaxis=dict(
                           title=ylabel,
                           titlefont=dict(
                               family='Courier New, monospace',
                               size=18,
                               color='#7f7f7f')))
        fig = go.Figure(data=data, layout=layout)
        return fig


    def test(self):
        a = [1,2,3,4,5]
        nv_min = 10
        nv_max = 20
        ov_min = min(a)
        ov_max = max(a)
        tmp_ = []

        scale = (nv_max - nv_min) / (ov_max - ov_min)
        for i in a:
            nw_ = (i - ov_min) * scale + nv_min
            tmp_.append(nw_)


    def get_scatter_plot_per_concentration(self, dataframe, xvalue, yvalue, spl_var, nnc, title, xlabel, ylabel,mut_label, *args):
        '''
        :param dataframe:
        :param xvalue:
        :param yvalue:
        :param std_err:
        :param spl_var:
        :param nnc:
        :param title:
        :param xlabel:
        :param ylabel:
        :param args:
        :return: plotly figure object
        '''
        # CONC_ = u'Assay Conc (nM)'
        # xvalue = u'BiAb AC-SINS (%)'
        # yvalue = "gmean"
        # std_err = "std"
        # nnc = "NNC"
        # mut_label = "FX Mut (HC|LC)"
        # title = "AC-SINS versus TGT normalized"
        # xlabel = 'AC-SINS score(%)'
        # ylabel = 'Ratio (WT)'
        data = []
        import types
        assert isinstance(args, types.TupleType)
        for i in args[0]:
            tmp_ = dataframe[(dataframe[spl_var] == i)]
            data.append(go.Scatter(x=tmp_[xvalue],
                                   y=tmp_[yvalue], \
                                   mode='markers', text=tmp_[nnc] + "\n" + \
                                                        tmp_[mut_label],
                                   name="conc: " + str(i) + " nM", marker=dict(size=10))
                        )
        layout = go.Layout(title=title,
                           xaxis=dict(title=xlabel,
                                      titlefont=dict(
                                          family='Courier New, monospace',
                                          size=18,
                                          color='#7f7f7f'
                                      )
                                      ),
                           yaxis=dict(
                               title=ylabel,
                               titlefont=dict(
                                   family='Courier New, monospace',
                                   size=18,
                                   color='#7f7f7f')))

        return go.Figure(data=data, layout=layout)



    def get_scatter_plot_w_error_bars_per_AA(self, dataframe, xvalue, yvalue, std_err, spl_var, nnc, title, xlabel, ylabel,
                                  mut_label, *args):
        data = []
        import types
        assert isinstance(args, types.TupleType)
        for i in args[0]:
            tmp_ = dataframe[(dataframe[spl_var].str.contains(i))]
            data.append(go.Scatter(x=tmp_[xvalue],
                               y=tmp_[yvalue], \
                               error_y=dict(type='data', array=tmp_[std_err], visible=True), \
                               mode='markers', text=tmp_[nnc] + "\n" + \
                                                    tmp_[mut_label],
                               name="conc: " + str(i) + " nM", marker=dict(size=10))
                    )
        layout = go.Layout(title=title,
                       xaxis=dict(title=xlabel,
                                  titlefont=dict(
                                      family='Courier New, monospace',
                                      size=18,
                                      color='#7f7f7f'
                                  )
                                  ),
                       yaxis=dict(
                           title=ylabel,
                           titlefont=dict(
                               family='Courier New, monospace',
                               size=18,
                               color='#7f7f7f')))
        return go.Figure(data=data, layout=layout)


    def get_scatter_plot_w_error_bars_m2(self, dataframe, yvalue, std_err, nnc, title, xlabel, ylabel,mut_label):
        data = []
        print dataframe.columns

        import types
        data.append(go.Scatter(y=dataframe[yvalue], \
                               error_y=dict(type='data', array=dataframe[std_err], visible=True), \
                               mode='markers', text=dataframe[nnc] + "\n" + \
                                                    dataframe[mut_label],marker=dict(size=10))
                    )
        layout = go.Layout(title=title,
                       xaxis=dict(title=xlabel,
                                  titlefont=dict(
                                      family='Courier New, monospace',
                                      size=18,
                                      color='#7f7f7f'
                                  )
                                  ),
                       yaxis=dict(
                           title=ylabel,
                           titlefont=dict(
                               family='Courier New, monospace',
                               size=18,
                               color='#7f7f7f')))
        return go.Figure(data=data, layout=layout)