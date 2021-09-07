#!/usr/bin/env python3

import graphviz
import pandas

# set paths
node_path = "node_list.xls"
edge_path = "edge_list.xls"

# settings
rm_external = True
rm_utility = True
rm_any = rm_external | rm_utility

# read data
node_df = pandas.read_excel(node_path)
edge_df = pandas.read_excel(edge_path)

# construct node mask
if rm_any:
    node_mask = False
    if rm_external:
        node_mask = node_mask | node_df['external_dependency'].values.astype("bool")
    if rm_utility:
        node_mask = node_mask | node_df['utility_dependency'].values.astype("bool")
    node_mask = ~node_mask
else:
    node_mask = slice(None)

# extract node data
node_submodule = node_df['submodule'].values[node_mask]
# node_module = node_df['module'].values[node_mask]
# node_clstype = node_df['smtype'].values[node_mask]
node_shape = node_df['shape'].values[node_mask]
node_style = node_df['style'].values[node_mask]
node_color = node_df['color'].values[node_mask]

# construct node mask
if rm_any:
    edge_mask = False
    if rm_external:
        edge_mask = edge_mask | edge_df['external_dependency'].values.astype("bool")
    if rm_utility:
        edge_mask = edge_mask | edge_df['utility_dependency'].values.astype("bool")
    edge_mask = ~edge_mask
else:
    edge_mask = slice(None)

# extract edge data
edge_submodule = edge_df['submodule'].values[edge_mask]
edge_dependency = edge_df['dependency'].values[edge_mask]
edge_label = edge_df['label'].values[edge_mask]

# create graph object
G = graphviz.Digraph('pybropt_arch', filename='pybropt_arch.gv')

# add nodes to graph
for name,shape,style,color in zip(node_submodule, node_shape, node_style, node_color):
    G.node(name, shape = shape, style = style, color = color)

# add edges to graph
for c,e in zip(edge_submodule, edge_dependency):
    G.edge(c, e)

G.view()
