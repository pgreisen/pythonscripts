from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from scipy.cluster.hierarchy import dendrogram, linkage
import pylab as plt
import pandas as pd
import igraph
import plotly.plotly as py
from plotly.graph_objs import *
import numpy as np

class PhyloTreeWActivity:


    def __init__(self):
        self.display_key = ""
        self.textlabelkey = "TextLabel"
        self.nv_min = 3
        self.nv_max = 6
        self.title = ""
        self.inner_node_color = '#ccccff'
        self.color_dict = {'RSK': '#e60000', 'SGK': '#ffff00', 'PKC': '#32cd32', 'DMPK': '#e600e6', \
                      'NDR': '#3366ff', 'GRK': '#8080ff', 'CSA': 'magenta', 'MAST': 'green', 'YANK': 'pink'}

    def set_intervals(self, df, label, column_id, mn, sd):
        '''
        interval of 3 with mn+/- sd
        :param df: dataframe
        :param label: the identifier for the compound
        :param df_column: column with values specified
        :return: lists with ids which has been split into 3 intervals
        '''
        # list with ids
        c1 = []
        c2 = []
        c3 = []

        ctff_1 = mn - sd
        ctff_3 = mn + sd

        for label_, val_ in zip(df[label], df[column_id]):
            if (val_ < ctff_1):
                c1.append(label_.strip())
            elif (val_ > ctff_3):
                c3.append(label_.strip())
            elif (val_ >= ctff_1 and val_ <= ctff_3):
                c2.append(label_.strip())
        return c1, c2, c3


    def convert_to_igraph(self,tree):
        # Convert a Biopython Tree object to an igraph Graph.
        def add_edge(graph, node1, node2):
            graph.add_edge(node1.name, node2.name)

        def build_subgraph(graph, top):
            """Traverse  the Tree, and retrieve  graph edges and nodes."""
            for clade in top:
                graph.add_vertex(name=clade.root.name)
                add_edge(graph, top.root, clade.root)
                build_subgraph(graph, clade)

        if tree.rooted:
            G = igraph.Graph(directed=True)
        else:
            G = igraph.Graph()
        G.add_vertex(name=str(tree.root))
        build_subgraph(G, tree.root)
        return G


    def get_tree(self,alignment, fileformat='fasta',distance='identity'):
        aln = AlignIO.read(alignment, fileformat)
        # Computations for the unrooted tree
        calculator = DistanceCalculator(distance)
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        return tree

    def set_node_color(self,label,node_color='blue'):
        for Key in self.color_dict:
            if Key in label:
                node_color = self.color_dict[Key]
        return node_color

    def set_node_color_w_experimental_data(self,V, dataframe, value1="NNC", value2="gmean"):
        '''

        :param V: vertex of graph
        :param df: dataframe
        :param label:
        :param label_id:
        :return:
        '''
        mn = dataframe[value2].mean()
        sd = dataframe[value2].std(ddof=1)
        # returns 3 lists with ids
        ctff_1, ctff_2, ctff_3 = self.set_intervals(dataframe, value1, value2, mn, sd)

        node_color = 'black'  # default color
        node_colors = []

        # loop over vertices
        for v in V:
            n_v_ = v['name']
            if ('Inner' in v['name']):
                node_colors.append("black")
            elif (n_v_ in ctff_1):
                node_colors.append("blue")
            elif (n_v_ in ctff_2):
                node_colors.append("white")
            elif (n_v_ in ctff_3):
                node_colors.append("red")
            else:
                node_color = 'yellow'#node_color  # default color
        return node_colors

    def get_display_labels(self, labels, dataframe, dataframe_key, *args):
        # format corerctly
        df_w_labels_ = {}
        # value2 = Seq ID
        display_labels_ = []
        index_ = dataframe[dataframe_key].values
        for label in labels:
            if 'Inner' in label:
                display_labels_.append('')
                continue
            else:
                pdbid_ = label
                if label in index_:
                    tmp_str = ""
                    for i in args[0]:
                        tmp_val_ = dataframe[dataframe[dataframe_key] == label][i].values[0]
                        if( isinstance(tmp_val_,float) ):
                            tmp_val_ = round(tmp_val_,2)
                        tmp_str = tmp_str + str(i) +": "+str(tmp_val_)+"\n"
                    tmp_label_ = str(dataframe_key)+" " + pdbid_ + "\n"+tmp_str
                    display_labels_.append(tmp_label_)
                else:
                    tmp_label_ = str(dataframe_key)+": " + pdbid_ + "\n NaN"
                    display_labels_.append(tmp_label_)
                df_w_labels_[label] = tmp_label_
        return display_labels_,pd.DataFrame(df_w_labels_.items(),columns=[dataframe_key,self.textlabelkey])


    def get_scaled_values(self,listofvalues):
        # print listofvalues
        ov_min = min(listofvalues)
        ov_max = max(listofvalues)
        tmp_ = []
        scale = (self.nv_max - self.nv_min) / (ov_max - ov_min)
        for i in listofvalues:
            nw_ = (i - ov_min) * scale + self.nv_min
            tmp_.append(nw_)
        return tmp_


    def debug(self,Xn, Yn, displaylabels):
        Xemp = {}
        Yemp = {}
        Xnode = {}
        Ynode = {}
        for x,y,z in zip(Xn,Yn,displaylabels):
            if(z != ''):
                # print x,y,z
                Xnode[z] = x
                Ynode[z] = y
            else:
                Xemp[z] = x
                Yemp[z] = y

        df_tmp_ = pd.merge(right=pd.DataFrame(Xemp.items(), columns=[self.display_key, 'Xemp']), right_on=self.display_key, \
                         left=pd.DataFrame(Yemp.items(), columns=[self.display_key, 'Yemp']), left_on=self.display_key)
        df_xy_ = pd.merge(right=pd.DataFrame(Xnode.items(), columns=[self.textlabelkey, 'Xnode']), right_on=self.textlabelkey, \
                         left=pd.DataFrame(Ynode.items(), columns=[self.textlabelkey, 'Ynode']), left_on=self.textlabelkey)
        return df_xy_


    def get_phylo_genetic_tree(self,sequencefile,dataframe,value1,value2,title, nodesizelabel, colorlabel,*args):
        self.title = title
        self.display_key = value1

        # alignment is fasta or clustalW format
        p_tree = self.get_tree(sequencefile)
        G = self.convert_to_igraph(p_tree)

        Edges = [e.tuple for e in G.es]
        V = [v for v in G.vs]

        node_colors = [self.inner_node_color if 'Inner' in v['name'] else self.set_node_color(v['name']) for v in V]

        labels = [v['name'] for v in V]

        display_labels = ['' if 'Inner' in label else label for label in labels]

        display_labels,display_w_labels = self.get_display_labels(labels, dataframe, value1, args[0])

        # merge datalabel with input data frame
        df_merge_ = pd.merge(right=display_w_labels,right_on=value1,left=dataframe,left_on=value1)

        layt = G.layout('kk')
        N = len(layt)  # N is equal to len(G.vs)

        Xn = [layt[k][0] for k in range(N)]
        Yn = [layt[k][1] for k in range(N)]

        Xe = []
        Ye = []

        for e in Edges:
            Xe += [layt[e[0]][0], layt[e[1]][0], None]
            Ye += [layt[e[0]][1], layt[e[1]][1], None]

        # merge dataframes for X,Y coordinates
        xynode_ = self.debug(Xn, Yn, display_labels)
        df_merge_ = pd.merge(right=display_w_labels, right_on=self.textlabelkey, left=xynode_, left_on=self.textlabelkey)
        df_merge_ = pd.merge(right=dataframe, right_on=value1,left=df_merge_,left_on=value1  )

        # edges
        trace1 = Scatter(x=Xe,
                         y=Ye,
                         mode='lines',
                         line=Line(color=self.inner_node_color, width=1),
                         hoverinfo='none'
                         )
        trace2 = Scatter(x=df_merge_["Xnode"],
                         y=df_merge_["Ynode"],
                         mode='markers',
                         name='',
                         marker=Marker(symbol='dot',
                                       size=df_merge_[nodesizelabel],
                                       sizemode='area',
                                       sizemin=4,
                                       sizeref=2*max(list(df_merge_[nodesizelabel]))/(35**2),
                                       colorscale='Viridis',
                                       showscale=True,
                                       colorbar=dict(title=colorlabel),
                                       color=df_merge_[colorlabel],
                                       line=Line(color='rgb(50,50,50)', width=0.5)
                                       ),

                         text=df_merge_["TextLabel"],
                         textposition='right',
                         textfont=dict(
                             family='sans serif',
                             size=8,
                             color='#ff7f0e'),
                         hoverinfo='text'
                         # hoverinfo='marker'
                         )
        '''
        trace3 = Scatter(x=xline_,
                         y=yline_,
                         mode='markers',
                         name='',
                         marker=Marker(symbol='dot',
                                       size=15, # node_size,
                                       line=Line(color='rgb(50,50,50)', width=0.5)
                                       ))
        '''
        axis = dict(showbackground=False,  #
                    showline=False,
                    zeroline=False,
                    showgrid=False,
                    showticklabels=False,
                    title=''
                    )
        layout = Layout(
            title=self.title,
            width=1000,
            height=1000,
            showlegend=True,
            xaxis=XAxis(axis),
            yaxis=YAxis(axis),
            margin=Margin(
                t=100
            ),
            hovermode='closest',
            annotations=Annotations([
                Annotation(
                    showarrow=False,
                    xref='paper',
                    yref='paper',
                    x=0.2,
                    y=-0.1,
                    xanchor='left',
                    yanchor='bottom',
                    font=Font(
                        size=14
                    )
                )
            ]),
        )

        data = Data([trace1, trace2]) # ,trace3])
        fig = Figure(data=data, layout=layout)
        return fig