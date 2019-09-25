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

import workflows
from workflows import settings as op2_settings
from workflows import colors

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
cseq3 = [blue, #petrol,
         turquoise, green, lime, orange, red, purple, lila]

cseq4 = [blue, turquoise, green, lime, yellow, orange, red, purple]

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
pa = PyPath()
tf = PyPath()
from pypath import data_formats
tf.init_network(data_formats.transcription)

with pypath.curl.cache_off():
    i = intercell.IntercellAnnotation()
    pa.init_network()

i.save_to_pickle(os.path.join(cachedir, 'intercell.pickle'))
pa.save_network(pfile=os.path.join(cachedir, 'network.pickle'))

#i = intercell.IntercellAnnotation(pickle_file=os.path.join(cachedir,
#                                                           'intercell.pickle'))

#pa.init_network(pfile=os.path.join(cachedir, 'network.pickle'))


print([x for x in dir(i) if not x.startswith('_')])
i.class_names
df = i.df
###############################################################################
## Checking ligand-receptor interactions
ligands = [p for p in i.classes['ligand'] if type(p) is str]
receptors = [p for p in i.classes['receptor'] if type(p) is str]

# Number of ligands and receptors
len(set(ligands))
len(set(receptors))

# Total connections
count = 0
lr_interacts = [pa.up_edge(lig, rec, directed=False) is not None
                for lig, rec in itertools.product(ligands, receptors)]
sum(lr_interacts)

r_per_l = []
for lig in ligands:
    aux = [pa.up_edge(lig, rec, directed=False) is not None for rec in receptors]
    r_per_l.append(sum(aux))

fig, ax = plt.subplots()

ax.hist(r_per_l, bins=100)
ax.set_title('Receptors per ligand')


l_per_r = []
for rec in ligands:
    aux = [pa.up_edge(lig, rec, directed=False) is not None for lig in ligands]
    l_per_r.append(sum(aux))

fig, ax = plt.subplots()

ax.hist(l_per_r, bins=100)
ax.set_title('Ligands per receptor')



## Checking number of ligands we have TF information for:
sum([lig in tf.vs['name'] for lig in ligands])

print('done')
###############################################################################

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
                                 index=[s.capitalize() for s in elem_by_class.\
                                 keys()]).sort_values()

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
[(k, len(set(v))) for k, v in elem_by_class.items()]
[(k, len(v)) for k, v in elem_by_class.items()]
# Connections between intercell classes

main_classes = list(elem_by_class.keys())
main_classes.remove('extracellular')
main_classes.remove('transmembrane')
main_classes.remove('secreted')
main_classes.remove('intracellular')
combs = list(itertools.combinations_with_replacement(main_classes, 2))
connections = dict((k, 0) for k in combs)
connections
for cat_a, cat_b in combs:
    print('Counting connections between %s and %s...' % (cat_a, cat_b))

    ups_a = [e for e in elem_by_class[cat_a] if not e.startswith('COMPLEX')]
    ups_b = [e for e in elem_by_class[cat_b] if not e.startswith('COMPLEX')]

    #print(ups_a[:3])
    for up_a, up_b in itertools.product(ups_a, ups_b):
        if pa.up_edge(up_a, up_b, directed=False):
            connections[(cat_a, cat_b)] += 1
#================================= PLOTTING ==================================#
# Proteins by class
fig, ax = plt.subplots()
#ax.barh(range(len(counts)), counts.values, color=blue)
ax.scatter(counts.values, range(len(counts)), color=blue)
ax.set_yticks(range(len(counts)))
labels = [s.replace('_', ' ').capitalize() for s in counts.index]
labels[labels.index('Ecm')] = 'ECM'
ax.set_yticklabels(labels)#, rotation=90)
ax.grid()
ax.set_title('Proteins by class')
ax.set_xscale('log')
ax.set_ylim(-1, len(counts))
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'intercell_prots_by_class.pdf'))

# Number of elements by class
fig, ax = plt.subplots()
rng = range(len(elem_counts_by_class))
ax.grid()
#ax.barh(rng, elem_counts_by_class.values, color=blue)
ax.scatter(elem_counts_by_class.values, rng, color=blue)
ax.set_yticks(rng)
labels = [s.replace('_', ' ').capitalize() for s in elem_counts_by_class.index]
labels[labels.index('Ecm')] = 'ECM'
ax.set_yticklabels(labels)
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

labels = [s.replace('_', ' ').capitalize() for s in elem_by_class.keys()]
labels[labels.index('Ecm')] = 'ECM'
cluster_hmap(sims, xlabels=labels, ylabels=labels,
             title='Szymkiewicz–Simpson similarity of major intercellular clas'
             + 'ses', filename=os.path.join(dest_dir,
                                            'intercell_similarity.pdf'))

# Entities by source
fig, ax = plt.subplots(figsize=(7, 7))
#ax.barh(range(len(elems_by_source)), elems_by_source.values, color=blue)
ax.grid()
ax.scatter(elems_by_source.values, range(len(elems_by_source)), color=blue)

ax.set_yticks(range(len(elems_by_source)))

labels = [s.capitalize() for s in elems_by_source.index]
labels[labels.index('Comppi')] = 'ComPPI'
labels[labels.index('Opm')] = 'OPM'
labels[labels.index('Hgnc')] = 'HGNC'
labels[labels.index('Hpmr')] = 'HPMR'
labels[labels.index('Topdb')] = 'TopDB'
labels[labels.index('Cellphonedb')] = 'CellPhoneDB'
labels[labels.index('Cspa')] = 'CSPA'
labels[labels.index('Matrixdb')] = 'MatrixDB'
labels[labels.index('Go')] = 'GO curated'

ax.set_yticklabels(labels)#, rotation=90)

ax.set_title('Entities by source')
#ax.set_xlim(-1, elems_by_source.max())
ax.set_ylim(-1, len(elems_by_source))
fig.tight_layout()
#ax.set_xscale('log')
fig.savefig(os.path.join(dest_dir, 'intercell_ents_by_source.pdf'))

# Some are empty!
i.classes['ligand_kirouac']

# Interactions between subclasses
edges = pd.DataFrame([[k[0], k[1], v] for (k, v) in connections.items()])

nodes = dict((k, len([x for x in elem_by_class[k] if not x.startswith('complex')]))
              for k in main_classes)
[(k, len(v)) for k, v in elem_by_class.items()]
nodes
nodes = pd.Series(nodes)
nodes.sort_index(inplace=True)
nodes
# "Adjacency"
adj = pd.DataFrame(index=nodes.index, columns=nodes.index)
for k, (s, t, v) in edges.iterrows():
    adj.loc[s, t] = v
    adj.loc[t, s] = v

labels = [s.replace('_', ' ').capitalize() for s in nodes.index]
labels[labels.index('Ecm')] = 'ECM'

#cmap = matplotlib.cm.get_cmap('jet')
#colors = list(map(cmap, np.linspace(1, 0, len(nodes))))

cmap = matplotlib.cm.get_cmap('gist_rainbow')
colors = list(map(cmap, np.linspace(0, 1, len(nodes))))
fig = chordplot(nodes, edges, #labels=labels, label_sizes=True,
                alpha=.5,
                colors=colors, figsize=(10, 8))#colors)
ax = fig.gca() # [left, bottom, width, height]
ax.set_position([0, 0, 0.75, 1]) # chordplot
ax2 = fig.add_axes([0.75, 0.75, 0.25, 0.25]) # legend
ax2.legend([matplotlib.lines.Line2D([0], [0], marker='o', color='w',
                                              markerfacecolor=c)
                      for c in colors],
                     labels,
                     loc='center', ncol=2)
ax2.set_axis_off()
ax3 = fig.add_axes([0.75, 0.5, 0.25, 0.25]) # barplot
ax4 =fig.add_axes([0.75, 0, 0.25, 0.5], sharex=ax3) # heatmap

ax3.bar(range(len(nodes)), nodes.values, color=colors)
#ax3.set_axis_off()
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.set_xticks([])
#ax.axes.get_yaxis().set_visible(False)

ax3.locator_params(axis='y', nbins=4)

#ax.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off', labelright='off', labelbottom='off')

#ax3.get_yaxis().set_visible(True)

mat = np.triu(adj.values).astype(float)
mat[np.where(mat == 0)] = np.nan


im = ax4.imshow(mat, cmap='plasma')
 # [left, bottom, width, height]
ax5 = fig.add_axes([0.75, 0.05, 0.25, 0.01]) # colorbar
fig.colorbar(im, cax=ax5, ax=ax5, orientation='horizontal')
ax4.set_axis_off()

fig
fig.savefig(os.path.join(dest_dir, 'intercell_interact_chordplot2.pdf'))

# =========================================================================== #
# Moving files to omnipath2_latex repository
tomove = [f for f in os.listdir(dest_dir)
          if (f.startswith('intercell') and f.endswith('.pdf'))]
# Remove the UpSet plot (too big)
tomove.remove('intercell_overlaps.pdf')

for f in tomove:
    shutil.copy2(os.path.join(dest_dir, f), os.path.join(latex_dir, f))
# =========================================================================== #
