#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os

import matplotlib.pyplot as plt
import pandas as pd

import pypath
from pypath import complex

# Colors!
cg = (87/255, 171/255, 39/255)
cb = (0/255, 84/255, 159/255)

cachedir = '/home/nico/pypath_cache'

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

pypath.settings.setup(cachedir=cachedir)
co = complex.ComplexAggregator()

# Number of complexes
total = len(co.complexes)
print('Currently present information for %d complexes' % total)


# Types of complex
homomultimer = 0
heteromultimer = 0
# Sources
unique_sources = set()
# References
unique_refs = set()

# Roll-it!
for k, v in co.complexes.items():
    unique_sources.update(v.sources)
    unique_refs.update(v.references)
    aux = k.split('-')
    unique_prots.update(aux)

    if len(aux) == 1:
        homomultimer += 1

    else:
        heteromultimer += 1

# Unique proteins across all complexes
print('Total number of unique proteins within complexes:', len(co.proteins))
print('Out of %d complexes, %d are homomultimers (%.2f %%) and %d '
      'heteromultimers (%.2f %%)' % (total, homomultimer,
                                     100 * homomultimer / total,
                                     heteromultimer,
                                     100 * heteromultimer / total))
fig, ax = plt.subplots()
ax.pie([homomultimer, heteromultimer],
       labels=['Homomultimers', 'Heteromultimers'],
       autopct='%.2f %%',
       colors=[cg, cb])
ax.set_title('Complex types')
fig.tight_layout()
fig.savefig('../figures/complex_types.svg')
fig
print('There are a total of %d unique references' % len(unique_refs))

# Complexes by resource
comp_by_res = dict(zip(unique_sources, [0] * len(unique_sources)))

# Roll-it again
for c in co.complexes.values():
    for s in c.sources:
        comp_by_res[s] += 1

comp_by_res = pd.Series(comp_by_res).sort_values(ascending=False)
rng = range(len(comp_by_res))
fig, ax = plt.subplots()
ax.bar(rng, comp_by_res.values, color=cb)
ax.set_xticks(rng)
ax.set_xticklabels(comp_by_res.index, rotation=90)
ax.set_xlabel('Resource')
ax.set_ylabel('Number of complexes')
fig.tight_layout()
fig.savefig('../figures/complex_by_source.svg')
