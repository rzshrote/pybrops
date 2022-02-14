#!/usr/bin/env python3

import graphviz
import pandas
import numpy

# set paths
graph_specs = "graph_specs.xls"

# settings
rm_external = True
rm_utility = True
rm_any = rm_external | rm_utility

# create graph object
G = graphviz.Digraph('pybrops_arch', filename='pybrops_arch.gv', format='png')
# G.attr(size='20,20!')
G.attr(dpi = "300")

################################################################################
################################################################################
################################################################################

# read node data
node_df = pandas.read_excel(graph_specs, sheet_name = "nodes")

# construct node mask
if rm_any:
    node_mask = False
    if rm_external:
        node_mask = node_mask | node_df['external_dependency'].values.astype("bool")
    if rm_utility:
        node_mask = node_mask | node_df['utility_dependency'].values.astype("bool")
    node_mask = node_mask | node_df['hide'].values.astype("bool")
    node_mask = ~node_mask
else:
    node_mask = slice(None)

# extract node data
node_submodule = node_df['submodule'].values[node_mask]
node_module = node_df['module'].values[node_mask]
# node_clstype = node_df['smtype'].values[node_mask]
node_shape = node_df['shape'].values[node_mask]
node_style = node_df['style'].values[node_mask]
node_color = node_df['color'].values[node_mask]

# get unique module names
node_module_uniq = numpy.unique(node_module)

# node clusters
node_cluster_context = {}
node_cluster_subgraph = {}
i = 0
for m in node_module_uniq:
    if m is not None:
        node_cluster_context[m] = G.subgraph(name = "cluster_{0}".format(i))
        node_cluster_subgraph[m] = node_cluster_context[m].__enter__()
        node_cluster_subgraph[m].attr(style='filled', color='lightgrey', label=m)
        i += 1

# add nodes to appropriate subgraph or graphs
for name,m,shape,style,color in zip(node_submodule, node_module, node_shape, node_style, node_color):
    if m is not None:
        node_cluster_subgraph[m].node(name, shape = shape, style = style, color = color)
    else:
        G.node(name, shape = shape, style = style, color = color)

# close contexts so that rendering is done properly
for k in node_cluster_context.keys():
    node_cluster_context[k].__exit__(None,None,None)

################################################################################
################################################################################
################################################################################

# read edge data
edge_df = pandas.read_excel(graph_specs, sheet_name = "edges")

# extract edge data
edge_submodule = edge_df['submodule'].values
edge_dependency = edge_df['dependency'].values
edge_label = edge_df['label'].values

# create a mask depending on what is in the graph
edge_mask1 = numpy.in1d(edge_submodule, node_submodule)
edge_mask2 = numpy.in1d(edge_dependency, node_submodule)
edge_mask = edge_mask1 & edge_mask2

# extract edges
edge_submodule = edge_submodule[edge_mask]
edge_dependency = edge_dependency[edge_mask]
edge_label = edge_label[edge_mask]

# add edges to graph
for c,e in zip(edge_submodule, edge_dependency):
    G.edge(c, e)

G.unflatten(stagger = 2)
G.view()
