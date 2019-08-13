#!/usr/bin/env python
#-*- coding: utf-8 -*-

#
#  This file is part of the `omnipath2` Python module
#
#  Copyright
#  2019
#  Heidelberg University, Uniklinik RWTH Aachen
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import cProfile as profile

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
proc.export_stats_by_category_pairs()
proc.export_annotations_by_entity()
proc.export_resources_by_entity()
proc.export_complexes_by_resource()

# testing memory usage of joins:

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

# profiling connection stats:

p = profile.Profile()
p.enable()
proc.export_stats_by_category_pairs()
p.disable()
