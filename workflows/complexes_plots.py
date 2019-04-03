#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os

import matplotlib.pyplot as plt
import pandas as pd

import pypath
from pypath import complex


#=================================== SETUP ===================================#

# Colors!
green = (87/255, 171/255, 39/255)
blue = (0/255, 84/255, 159/255)

# Setting up working environment
cachedir = '/home/nico/pypath_cache'

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

pypath.settings.setup(cachedir=cachedir)

#============================== RETRIEVING INFO ==============================#

co = complex.ComplexAggregator()

# Number of complexes
total = len(co.complexes)
print('Currently present information for %d complexes' % total)

# 1) Complex types
homomultimer = 0
heteromultimer = 0
# All sources
unique_sources = set()
# All references
unique_refs = set()

# Roll-it!
for k, v in co.complexes.items():
    unique_sources.update(v.sources)
    unique_refs.update(v.references)
    aux = k.split('-')

    if len(aux) == 1:
        homomultimer += 1

    else:
        heteromultimer += 1

homopt = 100 * homomultimer / total
heteropt = 100 * heteromultimer / total

# Unique proteins across all complexes
print('Total number of unique proteins within complexes:', len(co.proteins))
print('Out of %d complexes, %d are homomultimers (%.2f %%) and %d '
      'heteromultimers (%.2f %%)' % (total, homomultimer, homopt,
                                     heteromultimer, heteropt))

print('There are a total of %d unique references' % len(unique_refs))

# 2) Complexes by resource

comp_by_res = dict(zip(unique_sources, [0] * len(unique_sources)))

# Roll-it again
for c in co.complexes.values():
    for s in c.sources:
        comp_by_res[s] += 1

comp_by_res = pd.Series(comp_by_res).sort_values(ascending=True)
#================================= PLOTTING ==================================#

# 1) Complex types
fig, ax = plt.subplots()
ax.pie([homomultimer, heteromultimer],
       labels=['Homomultimers\n%.2f %%' % homopt,
               'Heteromultimers\n%.2f %%' % heteropt],
       colors=[green, blue])
ax.set_title('Complex types')
fig.tight_layout()
fig.savefig('../figures/complex_types.svg')

# 2) Complexes by resource
rng = range(len(comp_by_res))
fig, ax = plt.subplots()
ax.barh(rng, comp_by_res.values, color=blue)
ax.set_yticks(rng)
ax.set_yticklabels(comp_by_res.index)
ax.set_ylabel('Resource')
ax.set_xlabel('Number of complexes')
fig.tight_layout()
fig.savefig('../figures/complex_by_source.svg')
