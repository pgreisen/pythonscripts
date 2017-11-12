import sys,shutil,os,subprocess,re


class PlotTree:

    def __init__(self):
        self.setup_requirements()

        self.clustalfile = ""
        # http://igraph.org/python/doc/igraph.Graph-class.html
        self.graph_layout = 'kk'
        self.title = "Title"
        self.inner_node_color = '#ccccff'
        self.color_dict = {'RSK': '#e60000', 'SGK': '#ffff00', 'PKC': '#32cd32', 'DMPK': '#e600e6',\
                           'NDR': '#3366ff','GRK': '#8080ff', 'CSA': 'magenta', 'MAST': 'green', 'YANK': 'pink'}





    def setup_requirements(self):
        igraph_py = "!pip install python-igraph"
        subprocess.Popen(igraph_py, shell=True).wait()
        import_igraph = "import igraph"
        subprocess.Popen(import_igraph, shell=True).wait()
        bio_1 = "from Bio import Phylo"
        bio_2 = "from Bio.Phylo.TreeConstruction import DistanceCalculator"
        bio_3 = "from Bio.Phylo.TreeConstruction import DistanceTreeConstructor"
        bio_4 = "from Bio import AlignIO"

        subprocess.Popen(bio_1, shell=True).wait()
        subprocess.Popen(bio_2, shell=True).wait()
        subprocess.Popen(bio_3, shell=True).wait()
        subprocess.Popen(bio_4, shell=True).wait()



    def set_node_color(self, label):
        node_color = 'blue'#default color
        for Key in color_dict:
            if Key in label:
                node_color = color_dict[Key]
        return node_color



    def set_node_color(self,label):
        label = label.split('_', 1)[-1]
        label = label.replace('+','_')
        node_color = 'yellow'#default color
        if(label in low):
            node_color = "blue"
        elif(label in medium):
            node_color = "white"
        elif(label in high):
            node_color = "red"
        else:
            node_color = 'yellow'#default color
        return node_color


    def convert_to_igraph(self, tree):
        '''
        Convert a Biopython Tree object to an igraph Graph.
        '''
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

    def plotly_setup(self):
        plotly_1 = "import plotly.plotly as py"
        plotly_2 = "from plotly.graph_objs import *"
        subprocess.Popen(plotly_1, shell=True)
        subprocess.Popen(plotly_2, shell=True)

    def setup_graph(self):
        Xn = [layt[k][0] for k in range(N)]
        Yn = [layt[k][1] for k in range(N)]

        Xe = []
        Ye = []
        for e in Edges:
            Xe += [layt[e[0]][0], layt[e[1]][0], None]
            Ye += [layt[e[0]][1], layt[e[1]][1], None]

        # edges
        trace1 = Scatter(x=Xe,
                         y=Ye,
                         mode='lines',
                         line=Line(color=inner_node_color, width=1),
                         hoverinfo='none'
                         )
        # nodes
        trace2 = Scatter(x=Xn,
                         y=Yn,
                         mode='markers+text',
                         # mode='markers',
                         name='',
                         marker=Marker(symbol='dot',
                                       size=node_size,
                                       color=node_colors,
                                       line=Line(color='rgb(50,50,50)', width=0.5)
                                       ),
                         # remove or add labels here
                         text='',  # display_labels,
                         textposition='right',
                         textfont=dict(
                             family='sans serif',
                             size=8,
                             color='#ff7f0e'),
                         hoverinfo=display_labels  # 'text'
                         )
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
            # showlegend=False,
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
                    text="Data source: ",
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
        data = Data([trace1, trace2])
        fig = Figure(data=data, layout=layout)


    def plot_interactive_tree(self):
        # py.sign_in('empet', 'api-key')
        import plotly.offline as offline
        import plotly.graph_objs as go
        from plotly.offline import plot_mpl
        name_of_tree = obs + "_tree.html"
        offline.init_notebook_mode()
        offline.iplot(fig, filename='unrooted-tree')
        offline.plot(fig, image='png', filename=name_of_tree, image_filename=name_of_tree + '_sequence_tree',auto_open=False)

    # get TGT values
    def get_tree(self,clustalfile):
        '''
        :param clustalfile:
        :return:tree object
        '''
        aln = AlignIO.read(self.clustalfile, 'clustal')
        # Computations for the unrooted tree
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        return tree


    def setup_igraph_tree(self,alignmentfile):
        '''
        :param alignmentfile:
        :return:
        '''
        tree = self.get_tree(alignmentfile)
        # next setup the igraph tree
        G = convert_to_igraph(tree)
        Edges = [e.tuple for e in G.es]
        V = [v for v in G.vs]
        node_colors = [inner_node_color if 'Inner' in v['name'] else set_node_color(v['name']) for v in V]
        labels = [v['name'] for v in V]
        display_labels = ['' if 'Inner' in label else label for label in labels]
        # node size is set here:
        node_size = [25 if d_label == 'CSA' else 12 for d_label in display_labels]
        layt = G.layout( self.graph_layout )
        N = len(layt)  # N is equal to len(G.vs)




    def test(self):
        pass













G = convert_to_igraph(tree)
Edges = [e.tuple for e in G.es]
V = [v for v in G.vs]
node_colors = [inner_node_color if 'Inner' in v['name'] else set_node_color(v['name']) for v in V]

labels = [v['name'] for v in V]

display_labels = ['' if 'Inner' in label else label for label in labels]

# node size is set here:
node_size = [25 if d_label == 'CSA' else 12 for d_label in display_labels]

layt = G.layout('kk')
N = len(layt)  # N is equal to len(G.vs)



