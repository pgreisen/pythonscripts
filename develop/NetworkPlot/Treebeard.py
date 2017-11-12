import networkx, pylab
from networkx.drawing.nx_agraph import graphviz_layout
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

import plotly.plotly as py
from plotly.graph_objs import *

import networkx as nx


# What color to give to the edges?
e_color = '#ccccff'
# What colors to give to the nodes with similar labels?
color_scheme = {'RSK': '#e60000', 'SGK': '#ffff00', 'PKC': '#32cd32', 'DMPK': '#e600e6', 'NDR': '#3366ff',
                'GRK': '#8080ff', 'PKA': 'magenta', 'MAST': 'green', 'YANK': 'pink'}
# What sizes to give to the nodes with similar labels?
size_scheme = {'RSK': 200, 'SGK': 150, 'PKC': 350, 'DMPK': 400, 'NDR': 280, 'GRK': 370, 'PKA': 325, 'MAST': 40,
               'YANK': 200}


# Edit this to produce a custom label to color mapping
def label_colors(label):
    color_to_set = 'blue'
    for label_subname in color_scheme:
        if label_subname in label:
            color_to_set = color_scheme[label_subname]
    return color_to_set


# Edit this to produce a custom label to size mapping
def label_sizes(label):
    # Default size
    size_to_set = 20
    for label_subname in size_scheme:
        if label_subname in label:
            size_to_set = size_scheme[label_subname]
    return size_to_set


# Draw a tree whose alignment is stored in msa.phy
def draw_tree():
    # This loads the default kinase alignment that should be in the same directory as this script
    aln = AlignIO.read('agc.aln', 'clustal')
    # This will construct the unrooted tree.
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    G = Phylo.to_networkx(tree)


    node_sizes = []
    labels = {}
    node_colors = []
    for n in G:
        label = str(n)
        if 'Inner' in label:
            # These are the inner tree nodes -- leave them blank and with very small sizes.
            node_sizes.append(1)
            labels[n] = ''
            node_colors.append(e_color)
        else:
            # Size of the node depends on the labels!
            node_sizes.append(label_sizes(label))
            # Set colors depending on our color scheme and label names
            node_colors.append(label_colors(label))
            # set the label that will appear in each node
            labels[n] = label
    # Draw the tree given the info we provided!
    pos = graphviz_layout(G)
    #print "graphviz: ",pos
    pos1 = nx.spring_layout(G)
    #print "spring_layout(G): ", pos1
    pos2 = nx.spectral_layout(G)
    #print "spectral_layout(G): ",pos2

    pos3 = nx.fruchterman_reingold_layout(G)
    #print "spring_layout(G): ",pos3

    # assert 1 == 0


    # pos = nx.spectral_layout(G)
    networkx.draw(G, pos, edge_color=e_color, node_size=node_sizes, labels=labels, with_labels=True,node_color=node_colors)
    #pylab.show()

    # pos = nx.spectral_layout(G)
    #networkx.draw(G, pos1, edge_color=e_color, node_size=node_sizes, labels=labels, with_labels=True,
    #              node_color=node_colors)
    #pylab.show()

    # pos = nx.spectral_layout(G)
    #networkx.draw(G, pos2, edge_color=e_color, node_size=node_sizes, labels=labels, with_labels=True,
    #              node_color=node_colors)
    #pylab.show()

    # pos = nx.spectral_layout(G)
    #networkx.draw(G, pos3, edge_color=e_color, node_size=node_sizes, labels=labels, with_labels=True,
    #              node_color=node_colors)
    #pylab.show()

    # which graph is the default graph

    pylab.savefig('example.png')

    # convert to plotly
    pos = nx.spring_layout(G)

    dmin = 1
    ncenter = 0
    for n in pos:
        x, y = pos[n]
        d = (x - 0.5) ** 2 + (y - 0.5) ** 2
        if d < dmin:
            ncenter = n
            dmin = d
    # ?????
    ##p = nx.single_source_shortest_path_length(G, ncenter)

    edge_trace = Scatter(
        x=[],
        y=[],
        line=Line(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    for edge in G.edges():
        print "edge", edge
        # print "New new new ", edge
        print G[edge[0]], G.node[edge[0]]
        print G[edge[1]], G.node[edge[1]]
        # assert 1 == 0
        x0, y0 = G.node[edge[0]]['pos']
        x1, y1 = G.node[edge[1]]['pos']
        edge_trace['x'] += [x0, x1, None]
        edge_trace['y'] += [y0, y1, None]

    node_trace = Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers',
        hoverinfo='text',
        marker=Marker(
            showscale=True,
            # colorscale options
            # 'Greys' | 'Greens' | 'Bluered' | 'Hot' | 'Picnic' | 'Portland' |
            # Jet' | 'RdBu' | 'Blackbody' | 'Earth' | 'Electric' | 'YIOrRd' | 'YIGnBu'
            colorscale='YIGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line=dict(width=2)))

    for node in G.nodes():
        x, y = G.node[node]['pos']
        node_trace['x'].append(x)
        node_trace['y'].append(y)
    for node, adjacencies in enumerate(G.adjacency_list()):
        node_trace['marker']['color'].append(len(adjacencies))
        node_info = '# of connections: ' + str(len(adjacencies))
        node_trace['text'].append(node_info)






if __name__ == '__main__':
    draw_tree()
