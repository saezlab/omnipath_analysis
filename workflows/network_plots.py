#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os
import shutil
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wordcloud

import pypath
from pypath.main import PyPath
import data_tools
from data_tools.plots import cluster_hmap
from data_tools.plots import upset_wrap
from data_tools.plots import venn
from data_tools.models import PowerLaw

#=================================== SETUP ===================================#

# Colors!
green = (87/255, 171/255, 39/255)
lime = (189/255, 205/255, 0/255)
blue = (0/255, 84/255, 159/255)
blue75 = (64/255, 127/255, 183/255)
blue50 = (142/255, 186/255, 229/255)
yellow = (255/255, 237/255, 0/255)
orange = (246/255, 168/255, 0/255)
petrol = (0/255, 97/255, 101/255)
turquoise = (0/255, 152/255, 161/255)
red = (161/255, 16/255, 53/255)
bordeaux = (161/255, 16/255, 53/255)
purple = (97/255, 33/255, 88/255)
lila = (122/255, 111/255, 172/255)

# Setting up working environment
cachedir = '/home/nico/pypath_cache'
dest_dir = '../figures'
latex_dir = '../../omnipath2_latex/figures'

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

pypath.settings.setup(cachedir=cachedir)

cseq = [blue, petrol, turquoise, green, lime, yellow, orange, red, bordeaux,
        lila, purple]

def color(*args, **kwargs):
    idx = np.random.randint(0, len(cseq))
    aux = [val * 255 for val in cseq[idx]]
    return "rgb({:.0f}, {:.0f}, {:.0f})".format(*aux)

#============================== RETRIEVING INFO ==============================#
pa = PyPath()

#pa.init_network()
pa.init_network(pfile=os.path.join(cachedir, 'network.pickle'))
pa.get_directed()

#pa.save_network(pfile=os.path.join(cachedir, 'network.pickle'))

# Undirected graph
pa.graph.vcount()
pa.graph.ecount()

# Edges by source
edge_sources = dict()
refs_per_edge = list(map(len, pa.graph.es['references']))

for e in pa.graph.es:
    for k in e['sources']:
        if k in edge_sources.keys():
            edge_sources[k].add(e.index)
        else:
            edge_sources[k] = set([e.index])

edge_sources = dict((k, len(v)) for (k, v) in edge_sources.items())
edge_sources = pd.Series(edge_sources).sort_values(ascending=True)

# Nodes by source
node_sources = dict()
refs_per_node = list(map(len, pa.graph.vs['references']))

for v in pa.graph.vs:
    for k in v['sources']:
        if k in node_sources.keys():
            node_sources[k].add(v.index)
        else:
            node_sources[k] = set([v.index])

node_sources = dict((k, len(v)) for (k, v) in node_sources.items())
node_sources = pd.Series(node_sources).sort_values(ascending=True)

#upset_wrap(node_sources)

# Signs and directions:
cats = ['directed', 'undirected', 'signed', 'unsigned']
sd = dict((k, set()) for k in cats)

for e in pa.graph.es:
    d = e['dirs']

    if d.is_directed():
        sd['directed'].add(e.index)

    else:
        sd['undirected'].add(e.index)

    if d.has_sign():
        sd['signed'].add(e.index)

    else:
        sd['unsigned'].add(e.index)


# SignaLink pathways
slk_path = dict()

for pws in pa.graph.vs['slk_pathways']:

    for pw in pws:

        # Unifying names
        pw = pw.upper().split('(')[0].strip().split('/')[0]

        if (pw == '' or pw == 'IIP'):
            continue

        elif pw == 'HH':
            pw = 'HEDGEHOG'

        if pw in slk_path.keys():
            slk_path[pw] += 1

        else:
            slk_path[pw] = 1

slk_path = pd.Series(slk_path).sort_values(ascending=True)

pa.graph.vs.attributes()

# Directed graph
#pa.dgraph.vcount()
#pa.dgraph.ecount()

# Node degrees and frequency - dict: keys = degree, values = counts
degrees = Counter(pa.graph.vs.degree())

#================================= PLOTTING ==================================#
# Edge sources
fig, ax = plt.subplots(figsize=(9, 7))
rng = range(len(edge_sources))
ax.barh(rng, edge_sources.values, color=blue)
ax.set_yticks(rng)
ax.set_yticklabels(edge_sources.index)
ax.set_ylabel('Resource')
ax.set_xlabel('Number of edges')
ax.set_xscale('log')
ax.grid(True, axis='x')
ax.set_ylim(-1, len(edge_sources))

fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'network_edges_by_source.pdf'))

# Edge sources - as wordcloud
wc = wordcloud.WordCloud(background_color='white', color_func=color)
wc.generate_from_frequencies(dict((k, np.log10(v))for (k, v) in edge_sources.iteritems()))

fig, ax = plt.subplots()
ax.imshow(wc)
ax.set_axis_off()
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'network_cloud_edges_by_source.pdf'))


# Node sources
fig, ax = plt.subplots(figsize=(9, 7))
rng = range(len(node_sources))
ax.barh(rng, node_sources.values, color=blue)
ax.set_yticks(rng)
ax.set_yticklabels(node_sources.index)
ax.set_ylabel('Resource')
ax.set_xlabel('Number of nodes')
ax.set_xscale('log')
ax.grid(True, axis='x')
ax.set_ylim(-1, len(node_sources))

fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'network_nodes_by_source.pdf'))

# Node sources - as wordcloud
wc = wordcloud.WordCloud(background_color='white', color_func=color)
wc.generate_from_frequencies(dict((k, np.log10(v))for (k, v) in node_sources.iteritems()))

fig, ax = plt.subplots()
ax.imshow(wc)
ax.set_axis_off()
fig.tight_layout()pa.graph.vs.degree()
fig.savefig(os.path.join(dest_dir, 'network_cloud_nodes_by_source.pdf'))

# Distribution of references by edge/node
fig, ax = plt.subplots()
bins = np.histogram(np.hstack([refs_per_edge, refs_per_node]), bins=100)[1]
ax.hist(refs_per_edge, color=blue, alpha=0.5, bins=bins, label='References per edge')
ax.hist(refs_per_node, color=green, alpha=0.5, bins=bins, label='References per node')
ax.set_yscale('log')
ax.legend()
ax.set_xlabel('Number of references')
ax.set_ylabel('Number of nodes/edges')
fig.tight_layout()

fig.savefig(os.path.join(dest_dir, 'network_refs_per_node_edge.pdf'))

# Distribution of node degrees
fig, ax = plt.subplots()
ax.hist(pa.graph.vs.degree(), color=blue, bins=100)
ax.set_yscale('log')

ax.set_xlabel('Node degree')
ax.set_ylabel('Number of nodes')
fig.tight_layout()

fig.savefig(os.path.join(dest_dir, 'network_node_degree.pdf'))

# Signs and directions:
venn(list(sd.values()), labels=list(sd.keys()), title='Direction and signs of '
     + 'network edges', c=[blue, green, red, yellow],
     filename=os.path.join(dest_dir, 'network_dirs_signs.pdf'))

# SignaLink pathways
fig, ax = plt.subplots()
rng = range(len(slk_path))
ax.barh(rng, slk_path.values, color=blue)
ax.set_yticks(rng)
ax.set_yticklabels(slk_path.index)
ax.set_ylabel('SignaLink pathway')
ax.set_xlabel('Members in network')
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'network_slk_pathw.pdf'))

# SignaLink pathways - as wordcloud
wc = wordcloud.WordCloud(background_color='white', color_func=color)
wc.generate_from_frequencies(slk_path.to_dict())

fig, ax = plt.subplots()
ax.imshow(wc)
ax.set_axis_off()
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'network_cloud_slk_pathw.pdf'))

# Node degrees vs. frequency

model = PowerLaw(list(degrees.keys()), list(degrees.values()))
fig = model.plot()
ax = fig.gca()
fmt_a = ('%.2E' % model.a).split('E')

ax.text(10, 1e4, r'$y= %s\times 10^{%d} x^{%.2f}$' % (fmt_a[0], int(fmt_a[1]), np.round(model.k, 2)), size=15)
ax.set_xlabel('Node degree')
ax.set_ylabel('Frequency')

fig.savefig(os.path.join(dest_dir, 'network_powerlaw.pdf'))

# =========================================================================== #

# Moving files to omnipath2_latex repository
tomove = [f for f in os.listdir(dest_dir)
          if (f.startswith('network') and f.endswith('.pdf'))]

for f in tomove:
    shutil.copy2(os.path.join(dest_dir, f), os.path.join(latex_dir, f))
# =========================================================================== #
