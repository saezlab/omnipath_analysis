#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os
import shutil

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import pypath
from pypath.main import PyPath
from pypath import annot

# Colors!
cg = (87/255, 171/255, 39/255)
cb = (0/255, 84/255, 159/255)

cachedir = '/home/nico/pypath_cache'
dest_dir = '../figures'
latex_dir = '../../omnipath2_latex/figures'

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

#============================== RETRIEVING INFO ==============================#
pypath.settings.setup(cachedir=cachedir)

#with pypath.curl.cache_off():
#    pa.init_network()
#    a = annot.AnnotationTable(keep_annotators=True, create_dataframe=True)
#    a.load()

#a.save_to_pickle(os.path.join(cachedir, 'annot.pickle'))

a = annot.AnnotationTable(pickle_file=os.path.join(cachedir, 'annot.pickle'))


df = a.to_dataframe()

print([x for x in dir(a) if not x.startswith('_')])

df.shape

# Number of annotations per protein/complex
annots_per_prot = [sum(df.iloc[i, :].values) for i in range(len(df))]

# Number of proteins by resource
prots_by_res = dict()

for c in df.columns:
    source = c.split('__')[0]
    prots = set(df.index[df[c]])

    if source not in prots_by_res.keys():
        prots_by_res[source] = prots

    else:
        prots_by_res[source].update(prots)

prots_by_res = pd.Series(list(map(len, prots_by_res.values())),
                         index=prots_by_res.keys()).sort_values()

# Number of annotation subclasses per resource
all_keys = np.array([i[0] for i in a.cols.keys()])
unique_keys = list(set(all_keys))
len(all_keys)
subclasses = pd.Series([sum(all_keys == x) for x in unique_keys],
                        index=unique_keys).sort_values()
a.annots.keys()
#================================= PLOTTING ==================================#
# Annotation classes by resource
rng = range(len(subclasses))
fig, ax = plt.subplots(figsize=(9, 7))
ax.grid(True, axis='x')
ax.barh(rng, subclasses.values, color=cb)
ax.set_yticks(rng)
ax.set_yticklabels(subclasses.index)
ax.set_ylim(-1, len(subclasses))
ax.set_xlabel('Number of annotation classes')
ax.set_xscale('log')
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'annot_classes_by_source.pdf'))

# Annotations per protein/complex
fig, ax = plt.subplots()
ax.hist(annots_per_prot, bins=100, color=cb)
ax.set_title('Annotations per protein/complex')
ax.set_xlabel('Number of annotations')
ax.set_ylabel('Proteins/complexes')
#ax.set_yscale('log')
ax.set_xlim([-1, max(annots_per_prot)])
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'annot_per_prot.pdf'))

# Proteins/complexes by resource
fig, ax = plt.subplots(figsize=(7, 6))
ax.grid(True, axis='x')
rng = range(len(prots_by_res))
ax.barh(rng, prots_by_res.values, color=cb)
ax.set_title('Proteins/complexes by resource')
ax.set_ylabel('Resource')
ax.set_xlabel('Number of proteins/complexes')
ax.set_ylim(-1, len(prots_by_res))
ax.set_yticks(rng)
ax.set_yticklabels(prots_by_res.index)
ax.set_xscale('log')
fig.tight_layout()
fig.savefig(os.path.join(dest_dir, 'annot_prot_by_source.pdf'))

# =========================================================================== #
# Moving files to omnipath2_latex repository
tomove = [f for f in os.listdir(dest_dir)
          if (f.startswith('annot') and f.endswith('.pdf'))]

for f in tomove:
    shutil.copy2(os.path.join(dest_dir, f), os.path.join(latex_dir, f))

#================================ WORKAROUND =================================#
###############################################################################
#           vvvv       HERE ON USES INFO FROM WEBSERVICE       vvvv           #
###############################################################################

#from urllib.request import urlopen
#from urllib.request import Request as Request

#from data_tools.databases import to_df

#url = 'http://omnipathdb.org/annotations'
#params = ['organisms=9606']
#data = '&'.join(params)

#req = Request(url)#'?&'.join([url, data]))

#response = urlopen(req)
#page = response.read(99999999).decode('utf-8')

#df = to_df(page, header=True)
#df.head()
#df.shape
#print('Total number of individual annotations:', len(set(df.record_id)))
#print('Total number of individual proteins (UniProt):', len(set(df.uniprot)))
#print('Total number of individual proteins (GeneSymbol):', len(set(df.genesymbol)))
#print('Total number of individual sources:', len(set(df.source)))

#ups_by_source = dict()

#for n, subdf in df.groupby('source'):
#    print(n + '\n' + '=' * len(n))
#    annots = set(subdf.label)
    #print(str(annots) + '\n')
#    for a in annots:
#        print(a)
#        subsubdf = subdf.loc[subdf.label == a, :]
#        print('\t' + str(set(subsubdf.value)))
#    ups_by_source[n] = set(subdf.uniprot)
#    print('\n')

#prots_by_ref = pd.Series(list(map(len, ups_by_source.values())),
#                         index=ups_by_source.keys())
#prots_by_ref.sort_values(ascending=False, inplace=True)

#fig, ax = plt.subplots()
#rng = range(len(prots_by_ref))
#ax.bar(rng, prots_by_ref, color=cb)

#ax.set_xticks(rng)
#ax.set_xticklabels(prots_by_ref.index, rotation=90)
#ax.set_yscale('log')
#ax.set_ylabel('Number of proteins')
#ax.set_xlabel('Source')
#fig.tight_layout()
#fig.savefig('../figures/annot_prot_by_source.svg')
#prots_by_ref
