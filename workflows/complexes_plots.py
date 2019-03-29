#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os

import matplotlib.pyplot as plt

import pypath
from pypath import complex

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists('new_cache'):
    os.makedirs('new_cache')

pypath.settings.setup(cachedir='new_cache')
co = complex.ComplexAggregator()

# Number of complexes
total = len(co.complexes)
print('Currently present information for %d complexes' % total)

# Unique proteins across all complexes
unique_prots = set()
# Types of complex
homomultimer = 0
heteromultimer = 0

for c in co.complexes.keys():
    aux = c.split('-')
    unique_prots.update(aux)

    if len(aux) == 1:
        homomultimer += 1

    else:
        heteromultimer += 1

print('Total number of unique proteins within complexes:', len(unique_prots))
print('Out of %d complexes, %d are homomultimers (%.2f %%) and %d '
      'heteromultimers (%.2f %%)' % (total, homomultimer,
                                     100 * homomultimer / total,
                                     heteromultimer,
                                     100 * heteromultimer / total))
fig, ax = plt.subplots()
ax.pie([homomultimer, heteromultimer],
       labels=['Homomultimers', 'Heteromultimers'],
       autopct='%.2f %%')
ax.set_title('Complex types')
fig.savefig('../figures/complex_types.svg')

co.complexes['P01584-P27930-Q9NPH3'].sources
