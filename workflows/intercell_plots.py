#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os
import itertools
import shutil

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.tight_layout import get_renderer
import pandas as pd
import wordcloud

import pypath
from pypath import intercell
from pypath.main import PyPath
#reload(data_tools)
#reload(data_tools.plots)
import data_tools
from data_tools.iterables import subsets
from data_tools.spatial import equidist_polar
from data_tools.iterables import similarity
from data_tools.plots import cluster_hmap
from data_tools.plots import upset_wrap
from data_tools.plots import chordplot

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

# Color sequences
cseq = [blue, blue75, petrol, turquoise, green, lime] # More gradual
cseq2 = [blue, green, orange, red, purple] # More contrasted
cseq3 = [blue, petrol, turquoise, green, yellow, orange, red, purple]

# Setting up the working environment
cachedir = '/home/nico/pypath_cache'
dest_dir = '../figures'
latex_dir = '../../omnipath2_latex/figures'

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

pypath.settings.setup(cachedir=cachedir)

#============================== RETRIEVING INFO ==============================#
#with pypath.curl.cache_off():
#    i = intercell.IntercellAnnotation()

#i.save_to_pickle(os.path.join(cachedir, 'intercell.pickle'))

i = intercell.IntercellAnnotation(pickle_file=os.path.join(cachedir,
                                                           'intercell.pickle'))

pa = PyPath()
#pa.init_network()

#pa.save_network(pfile=os.path.join(cachedir, 'network.pickle'))
pa.init_network(pfile=os.path.join(cachedir, 'network.pickle'))



print([x for x in dir(i) if not x.startswith('_')])
i.class_names
df = i.df

# Elements by class
elem_by_class = dict()

for c in set(df.mainclass):
    elems = set(df.loc[df.mainclass == c, 'uniprot'])

    if pd.isna(c):
        continue

    if c not in elem_by_class.keys():
        elem_by_class[c] = elems

    else:
        elem_by_class[c].update(elems)

elem_counts_by_class = pd.Series(map(len, elem_by_class.values()),
                                 index=elem_by_class.keys()).sort_values()

# Number of proteins by class
counts = dict((c, i.counts()[c]) for c in i.class_names)
counts = pd.Series(counts).sort_values(ascending=True)

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

# Connections between intercell classes

aux = list(elem_by_class.keys())
aux.remove('extracellular')
aux.remove('transmembrane')
aux.remove('secreted')

combs = list(itertools.combinations_with_replacement(aux, 2))
connections = dict((k, 0) for k in combs)

for cat_a, cat_b in combs:
    print('Counting connections between %s and %s...' % (cat_a, cat_b))

    ups_a = [e for e in elem_by_class[cat_a] if not e.startswith('COMPLEX')]
    ups_b = [e for e in elem_by_class[cat_b] if not e.startswith('COMPLEX')]

    for up_a, up_b in itertools.product(ups_a, ups_b):
        if pa.up_edge(up_a, up_b):
            connections[(cat_a, cat_b)] += 1

#================================= PLOTTING ==================================#
# Proteins by class
fig, ax = plt.subplots()
ax.barh(range(len(counts)), counts.values, color=blue)
ax.set_yticks(range(len(counts)))
ax.set_yticklabels([s.replace('_', ' ') for s in counts.index])#, rotation=90)
ax.set_title('Proteins by class')
ax.set_xscale('log')
ax.set_ylim(-1, len(counts))
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'intercell_prots_by_class.pdf'))

# Number of elements by class
fig, ax = plt.subplots()
rng = range(len(elem_counts_by_class))
ax.grid(True, axis='x')
ax.barh(rng, elem_counts_by_class.values, color=blue)
ax.set_yticks(rng)
ax.set_yticklabels([s.replace('_', ' ') for s in elem_counts_by_class.index])
ax.set_title('Number of elements per intercell class')
ax.set_xscale('log')
ax.set_ylim(-1, len(elem_counts_by_class))
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'intercell_prots_by_class2.pdf'))

# By similarity
groups = list(elem_by_class.values())
sims = []

for (a, b) in itertools.product(groups, repeat=2):
    sims.append(similarity(set(a), set(b), mode='ss'))

sims = np.array(sims).reshape(len(groups), len(groups))

labels = [s.replace('_', ' ') for s in elem_by_class.keys()]

cluster_hmap(sims, xlabels=labels, ylabels=labels,
             title='Szymkiewicz–Simpson similarity of major intercellular clas'
             + 'ses', filename=os.path.join(dest_dir,
                                            'intercell_similarity.pdf'))

# Entities by source
fig, ax = plt.subplots(figsize=(7, 7))
ax.barh(range(len(elems_by_source)), elems_by_source.values, color=blue)
ax.set_yticks(range(len(elems_by_source)))
ax.set_yticklabels(elems_by_source.index)#, rotation=90)
ax.set_title('Entities by source')
ax.set_xscale('log')
ax.set_ylim(-1, len(elems_by_source))
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'intercell_ents_by_source.pdf'))

# Some are empty!
i.classes['ligand_kirouac']

# Interactions between subclasses
edges = pd.DataFrame([[k[0], k[1], v] for (k, v) in connections.items()])
nodes = dict((k, len([x for x in elem_by_class[k] if not x.startswith('complex')]))
              for k in aux)
nodes = pd.Series(nodes)
nodes.sort_index(inplace=True)

labels = [s.replace('_', ' ') for s in nodes.index]

#cmap = matplotlib.cm.get_cmap('jet')
#colors = list(map(cmap, np.linspace(1, 0, len(nodes))))
fig = chordplot(nodes, edges, labels=labels, label_sizes=True, alpha=0.5, colors=cseq3)#colors)
fig.savefig(os.path.join(dest_dir, 'intercell_interact_chordplot.pdf'))
# =========================================================================== #
# Moving files to omnipath2_latex repository
tomove = [f for f in os.listdir(dest_dir)
          if (f.startswith('intercell') and f.endswith('.pdf'))]
# Remove the UpSet plot (too big)
tomove.remove('intercell_overlaps.pdf')

for f in tomove:
    shutil.copy2(os.path.join(dest_dir, f), os.path.join(latex_dir, f))
# =========================================================================== #
