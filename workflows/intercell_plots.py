#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os

import matplotlib.pyplot as plt
import pandas as pd

import pypath
from pypath import intercell
from pypath.main import PyPath
from data_tools.plots import venn

# Colors!
cg = (87/255, 171/255, 39/255)
cb = (0/255, 84/255, 159/255)

# TODO: Vertical barplots
# TODO: Colors!!

cachedir = '/home/nico/pypath_cache'

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

pypath.settings.setup(cachedir=cachedir)

i = intercell.IntercellAnnotation()

# Proteins by resource
elems = [a for a in dir(i) if a.endswith('_by_resource')]
resources = set()

for e in elems:
    resources.update(getattr(i, e).keys())

resources = list(resources)

df = pd.DataFrame(0, columns=resources, index=elems)

for e in elems:
    aux = getattr(i, e)
    for k, v in aux.items():
        df.loc[e, k] = len(v)

df.sum(axis=1)

fig, ax = plt.subplots()
rng = range(len(df))

df = df.loc[['receptors_by_resource', 'ecm_by_resource', 'ligands_by_resource']]
df
for (idx, typ) in enumerate(df.columns):

    if idx == 0:
        btm = df.loc[:, typ].copy().values
        ax.bar(rng, btm, label=typ)

    else:
        tmp = df.loc[:, typ].copy().values
        ax.bar(rng, tmp, label=typ, bottom=btm)
        btm += tmp
df
ax.set_xticks(rng)
ax.set_xticklabels(['%s (%d)' %(k.split('_')[0], df.loc[k, :].sum()) for k in df.index], rotation=90)
ax.legend()
fig.tight_layout()
fig.savefig('../figures/intercel_by_source.svg')

# Overlaps by source
for e in elems:
    aux = getattr(i, e)
    plot = venn(list(aux.values()), list(aux.keys()), sizes=True)
    ax = plot.gca()
    name = e.split('_')[0]
    ax.set_title('Overlap of %s by resource' % name)
    plot.savefig('../figures/intercel_%s_by_source.svg' % name)

# Interactions and overlaps between intercellular annotations
pa = PyPath()
pa.init_network()


elems = ['ecm_proteins', 'receptors', 'ligands']

annots = dict(zip(elems, [getattr(i, e) for e in elems]))
venn(list(annots.values()), list(annots.keys()), sizes=True, )

##pa.get_edge('P01133', 'P00533')
