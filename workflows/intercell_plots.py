#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os
import itertools as itt

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

import pypath
from pypath import intercell
from pypath.main import PyPath
from data_tools.plots import venn, upset_wrap
from data_tools.spatial import equidist_polar


#=================================== SETUP ===================================#
# Colors!
green = (87/255, 171/255, 39/255)
lime = (189/255, 205/255, 0/255)
blue = (0/255, 84/255, 159/255)
blue75 = (64/255, 127/255, 183/255)
orange = (246/255, 168/255, 0/255)
petrol = (0/255, 97/255, 101/255)
turquoise = (0/255, 152/255, 161/255)
red = (161/255, 16/255, 53/255)
purple = (97/255, 33/255, 88/255)

# Color sequences
cseq = [blue, blue75, petrol, turquoise, green, lime] # More gradual
cseq2 = [blue, green, orange, red, purple] # More contrasted

# Setting up the working environment
cachedir = '/home/nico/pypath_cache'

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

pypath.settings.setup(cachedir=cachedir)

#============================== RETRIEVING INFO ==============================#
with pypath.curl.cache_off():
    i = intercell.IntercellAnnotation()

print([x for x in dir(i) if not x.startswith('_')])

# Number of proteins by class
counts = dict((c, i.counts()[c]) for c in i.class_names)
counts = pd.Series(counts).sort_values(ascending=True)

# Overlapping proteins across classes


# Entities by source
sources = set([x.split('_')[-1] for x in i.classes.keys() if
(x not in i.class_names and not len(x.split('_'))==1)])

elems_by_source = dict()

for k, v in i.classes.items():
    s = k.split('_')[-1]

    if s in sources:

        if s not in elems_by_source.keys():
            elems_by_source[s] = v

        else:
            elems_by_source[s].update(v)

elems_by_source = pd.Series(map(len, elems_by_source.values()),
                            index=elems_by_source.keys())
elems_by_source.sort_values(inplace=True)

#================================= PLOTTING ==================================#
# Proteins by class
fig, ax = plt.subplots()
ax.barh(range(len(counts)), counts.values, color=blue)
ax.set_yticks(range(len(counts)))
ax.set_yticklabels(counts.index)#, rotation=90)
ax.set_title('Proteins by class')
fig.tight_layout()
fig.savefig('../figures/intercell_prots_by_class.svg')

# Entities by source
fig, ax = plt.subplots(figsize=(7, 7))
ax.barh(range(len(elems_by_source)), elems_by_source.values, color=blue)
ax.set_yticks(range(len(elems_by_source)))
ax.set_yticklabels(elems_by_source.index)#, rotation=90)
ax.set_title('Entities by source')
fig.tight_layout()
fig.savefig('../figures/intercell_ents_by_source.svg')

###############################################################################
#           vvvv       HERE ON USED FOR PREVIOUS VERSION       vvvv           #
###############################################################################

# List of intercellular attributes (ordered by resource)
#elems_s = [a for a in dir(i) if a.endswith('_by_resource')]

# List of unique resources for intercellular information
#resources = set()

#for e in elems_s:
#    resources.update(getattr(i, e).keys())

#resources = list(resources)

# 1) Number of proteins per resource and type (e.g. ecm, receptor, ligand)
#df = pd.DataFrame(0, columns=resources, index=elems_s)

#for e in elems_s:
#    aux = getattr(i, e)
#    for k, v in aux.items():
#        df.loc[e, k] = len(v)

# Reordering
#df = df.loc[['ligands_by_resource',
#             'ecm_by_resource',
#             'receptors_by_resource']]

# 3) Overlaps between intercellular annotations
#dir(i)

#elems = ['ecm_proteins', 'receptors', 'ligands']

annots = i.classes#dict(zip(elems, [getattr(i, e) for e in elems]))

# 4) Interactions between annotation types
with pypath.curl.cache_off():
    pa = PyPath()
    pa.init_network(redownload=True)

interact = dict()

# For all possible pairs of annotation types
for (a, b) in itt.combinations(annots.keys(), 2):
    # Check for edges in PPI network for any combination of proteins
    edges = [(1 if pa.get_edge(i, j) else 0)
             for (i, j) in itt.product(annots[a], annots[b])]
    # Store number of connecting edges
    interact[(a, b)] = sum(edges)

# Make a graph
G = nx.Graph()

for (i, (k, v)) in enumerate(annots.items()):
    G.add_node(k, size=len(v), color=cseq2[i])

for ((a, b), v) in interact.items():
    G.add_edge(a, b, weight=v)

#================================= PLOTTING ==================================#

# 1) Cummulative barplot of proteins per type and resource
#fig, ax = plt.subplots()
#rng = range(len(df))

#for (idx, typ) in enumerate(df.columns):

#    if idx == 0:
#        btm = df.loc[:, typ].copy().values
#        ax.barh(rng, btm, label=typ, color=cseq[idx])

#    else:
#        tmp = df.loc[:, typ].copy().values
#        ax.barh(rng, tmp, label=typ, left=btm, color=cseq[idx])
#        btm += tmp

#ax.set_yticks(rng)
#ax.set_yticklabels(['%s (%d)' %(k.split('_')[0], df.loc[k, :].sum())
#                    for k in df.index])
#ax.legend()
#fig.tight_layout()
#fig.savefig('../figures/intercel_by_source.svg')

# 2) Overlaps by source
#for e in elems_s:
#    aux = getattr(i, e)
#    plot = venn(list(aux.values()), list(aux.keys()), sizes=True, c=cseq2)
#    ax = plot.gca()
#    name = e.split('_')[0]
#    ax.set_title('Overlap of %s by resource' % name)
#    plot.savefig('../figures/intercel_%s_by_source.svg' % name)

# 3) Overlapping across annotations
#plot = venn(list(annots.values()), list(annots.keys()), sizes=True,
#            c=cseq2, filename='../figures/intercell_annots.svg')

# 4) Interactions between annotation types
nsizes = [G.nodes[n]['size'] * 2 for n in G.nodes]
ncolor = [G.nodes[n]['color'] for n in G.nodes]
esizes = [G.edges[e]['weight'] / 100 for e in G.edges]

fig, ax = plt.subplots(figsize=(7, 7))

# Get equidistant positions
pos = dict(zip(annots.keys(), equidist_polar(len(annots), r=0.25)))

nx.draw_networkx_edges(G, pos, width=esizes, ax=ax, alpha=0.5)
nx.draw_networkx_nodes(G, pos, node_color=ncolor, node_size=nsizes, ax=ax)
nx.draw_networkx_labels(G, pos, ax=ax)
nx.draw_networkx_edge_labels(G, pos, ax=ax,
                             edge_labels=dict(zip(G.edges,
                                                  [int(100 * e)
                                                   for e in esizes])))
ax.set_axis_off()
fig.tight_layout()

fig.savefig('../figures/intercel_edges.svg')
