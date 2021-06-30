#!/usr/bin/env python3
import pandas
import igraph
# import matplotlib.pyplot as plt

# read data
df = pandas.read_csv("architecture_graph.csv")

# generate edge list
ls = list(zip(df["class"],df["extends"]))

# create graph
G = igraph.Graph.TupleList(ls, directed = True)

# create visual style
visual_style = {}
visual_style["vertex_size"] = 50
visual_style["vertex_label"] = G.vs["name"]
visual_style["vertex_shape"] = "rectangle"
visual_style["label_size"] = 8
# visual_style["edge_curved"] = True
# visual_style["layout"] = G.layout("kamada_kawai")
visual_style["layout"] = G.layout("fruchterman_reingold")
# visual_style["layout"] = G.layout("reingold_tilford")
# visual_style["layout"] = G.layout("circle")
visual_style["bbox"] = (1000, 1000)
visual_style["margin"] = 100

# create plot as pdf document
igraph.plot(G, "architecture.pdf", **visual_style)
