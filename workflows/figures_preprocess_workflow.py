#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Denes Turei 2019
# turei.denes@gmail.com

import numpy as np
import pandas as pd

from workflows import figures_preprocess as figpreproc

omnipath_pickle = 'omnipath.pickle'
annotation_pickle = 'annotations2.pickle'
complex_pickle = 'complexes.pickle'
intercell_pickle = 'intercell.pickle'

proc = figpreproc.FiguresPreprocess(
    omnipath_pickle = omnipath_pickle,
    annotation_pickle = annotation_pickle,
    complex_pickle = complex_pickle,
)

proc.setup()
proc.load()

n = proc.network.records
i = proc.intercell.df

n_i = np.random.randint(0, high = n.shape[0] - 1, size = 100000)
i_i = np.random.randint(0, high = i.shape[0] - 1, size = 190000)

#n_s = n.loc[n_i,:]
n_s = n
# i_s = i.loc[i_i,:]
i_s = i

i_n = pd.merge(
    n_s,
    i_s,
    suffixes = ['', '_a'],
    how = 'inner',
    left_on = 'id_a',
    right_on = 'uniprot',
)
i_n.id_a = i_n.id_a.astype('category')

i_n_i = pd.merge(
    i_n,
    i_s,
    suffixes = ['_a', '_b'],
    how = 'inner',
    left_on = 'id_b',
    right_on = 'uniprot',
)
i_n_i.id_b = i_n_i.id_b.astype('category')

i_n_i.info(memory_usage = 'deep')

