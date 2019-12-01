import pygraphviz as pgv

G = pgv.AGraph()


G.add_node('a')
G.add_edge('b','c')
G.add_edge('b','d')
G.layout()


G.draw('file.png')
