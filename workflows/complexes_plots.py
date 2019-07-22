#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os
import shutil

import matplotlib.pyplot as plt
import pandas as pd

import pypath
from pypath import complex
from data_tools.plots import upset_wrap

#=================================== SETUP ===================================#

# Colors!
green = (87/255, 171/255, 39/255)
blue = (0/255, 84/255, 159/255)

# Setting up working environment
cachedir = '/home/nico/pypath_cache'
dest_dir = '../figures'
latex_dir = '../../omnipath2_latex/figures'

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

pypath.settings.setup(cachedir=cachedir)

#============================== RETRIEVING INFO ==============================#
with pypath.curl.cache_off():
    co = complex.ComplexAggregator()#resources=['Signor', 'Corum', 'CellPhoneDB',
                                    #          'Havugimana', 'Compleat',
                                    #          'ComplexPortal', 'Pdb', 'Hpmr',
                                    #          'GuideToPharmacology'])
                                              #, 'Humap'])

co.make_df()
sum(co.df.duplicated())

[i for i in dir(co) if not i.startswith('_')]

dir(co.complexes[list(co.complexes.keys())[0]])
len(co.complexes[list(co.complexes.keys())[0]].components)

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
    #aux = k.split('-')

    if len(v.components) == 1:
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

# 3) Overlaps across resources

len(co.df.components) == len(set(co.df.components))

comp_in_res = dict()
row_srcs = [r.split(';') for r in co.df.sources]

for i, srcs in enumerate(row_srcs):
    for s in srcs:
        if s not in comp_in_res.keys():
            comp_in_res[s] = set([i])
        else:
            comp_in_res[s].add(i)

#================================= PLOTTING ==================================#

# 1) Complex types
fig, ax = plt.subplots()
ax.bar([0, 1], [homomultimer, heteromultimer], color=blue)
ax.set_title('Complex types')
ax.set_xticks([0, 1])
ax.set_xticklabels(['Homomultimers\n%.2f %%' % homopt,
                    'Heteromultimers\n%.2f %%' % heteropt])
ax.text(0, 2000, '%d' % homomultimer, color='w', ha='center')
ax.text(1, 9000, '%d' % heteromultimer, color='w', ha='center')
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'complex_types.pdf'))

# 2) Complexes by resource
rng = range(len(comp_by_res))
fig, ax = plt.subplots()
ax.barh(rng, comp_by_res.values, color=blue)
ax.set_yticks(rng)
ax.set_yticklabels(comp_by_res.index)
ax.set_ylabel('Resource')
ax.set_xlabel('Number of complexes')
ax.set_xscale('log')
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'complex_by_source.pdf'))

# 3) Overlaps across resources
help(upset_wrap)
#upset_wrap(list(comp_in_res.values()), labels=comp_in_res.keys(),
#           drop_empty=True)
# =========================================================================== #
# Moving files to omnipath2_latex repository
tomove = [f for f in os.listdir(dest_dir)
          if (f.startswith('complex') and f.endswith('.pdf'))]

for f in tomove:
    shutil.copy2(os.path.join(dest_dir, f), os.path.join(latex_dir, f))
# =========================================================================== #
