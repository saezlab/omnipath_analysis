#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019 Saez Lab
#
# OmniPath2 analysis and figures suit
#
# Authors:
#
# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de
#
# Dénes Türei
# turei.denes@gmail.com
#

import collections

from pypath import settings
from pypath import common


_defaults = {
    
    # pickle dumps of all databases
    'pickle_dir': 'pickles',
    
    # directory for exported tables
    'tables_dir': 'tables',
    
    # directory for figures
    'figures_dir': 'figures',
    
    # directory with palettes
    'palette_dir': 'palettes',
    
    'timestamp_dirs': True,
    
    # pickles
    'omnipath_pickle': 'omnipath_pw_es_lr_20190918.pickle',
    'curated_pickle': 'network_curated_20190924.pickle',
    'complex_pickle': 'complexes.pickle',
    'annotations_pickle': 'annotations2.pickle',
    'intercell_pickle': 'intercell.pickle',
    'enz_sub_pickle': 'ptms.pickle',
    'tf_target_pickle': 'tfregulons-%s.pickle',
    'tf_mirna_pickle': 'tfmirna.pickle',
    'mirna_mrna_pickle': 'mirna_mrna.pickle',

    # supplementary tables
    'omnipath_tsv': 'network_summaries.tsv',
    'curated_tsv': 'network-curated_summaries.tsv',
    'tf_target_tsv': 'network-tf.tsv',
    'mirna_mrna_tsv': 'network-mirna-target_summaries.tsv',
    'tf_mirna_tsv': 'network-tf-mirna_summaries.tsv',
    'annotations_tsv': 'annotations.tsv',
    'enz_sub_tsv': 'enzsub_summaries.tsv',
    'intercell_tsv': 'intercell_summaries.tsv',
    'complex_tsv': 'complex_summaries.tsv',
    
    # for R plotting
    'intercell_classes_tsv': 'intercell_classes_%s.tsv',
    'main_coverage_tsv': 'main_coverage_%s.tsv',
    'intercell_network_by_resource_tsv': 'stats_by_resource_%s.tsv',
    'resources_by_entity_tsv': 'resources_by_entity_%s.tsv',
    'complexes_by_resource_tsv': 'complexes_by_resource_%s.tsv',
    'connections_tsv': 'connections_%s.tsv',
    'category_overlaps_tsv': 'category_overlaps_%s.tsv',

    # tfregulons levels
    'tfregulons_levels': {'A', 'B', 'C', 'D'},

    # datasets
    'datasets': [
       'omnipath',
       'curated',
       'complex',
       'annotations',
       'intercell',
       'tf_target',
       'tf_mirna',
       'mirna_mrna',
       'enz_sub',
    ],
    
    'omnipath_mod': 'main',
    'curated_mod': 'main',
    'complex_mod': 'complex',
    'annotations_mod': 'annot',
    'intercell_mod': 'intercell',
    'enz_sub_mod': 'ptm',
    'tf_target_mod': 'main',
    'tf_mirna_mod': 'main',
    'mirna_mrna_mod': 'main',
    
    'omnipath_args': {
        'use_omnipath': True,
        'kinase_substrate_extra': True,
        'ligand_receptor_extra': True,
        'pathway_extra': True,
    },
    
    'dependencies': {
        'intercell': ('annotations',),
        'annotations': ('complex',),
    },
    
    # figures graphic param defaults
    'font_family': [
        # we want to use this:
        'Helvetica Neue LT Std',
        # this is only fallback:
        'sans',
    ],
    'font_variant': 'normal',
    'font_style': 'normal',
    'font_stretch': 'condensed',
    'font_weight': 'medium',
    'font_size': 12,
    'font_sizes': {
        'axis_label': 1.,
        'ticklabel': .8,
        'legend_title': 1.,
        'legend_label': .8,
        'title': 1.2,
        'annotation': .8,
    },
    'fig_height': 3,
    'fig_width': 4,
    
    # figure filenames
    'inter_class_degree_pdf': 'inter_class_degree_%s_%s_%s_%s',
    'counts_by_class_pdf': 'counts_by_class_%s_%s',
    'counts_by_resource_pdf': 'counts_by_resource_%s',
    'inter_class_sim_pdf': 'inter_class_sim',
    'inter_class_chordplot_pdf': 'inter_class_chordplot_%s_%s_%s',
    
    # default palette
    'palette': 'Nico_3',
    
}


nico = {
    'environment': 'nico',
    'cachedir': '/home/nico/pypath_cache',
    'figures_dir': 'figures',
    'latex_dir': '../omnipath2_latex/figures',
}

denes = {
    'environment': 'denes',
    'latex_dir': 'latex',
    'figures_dir': 'figures',
}



def reset_all():
    
    settings = collections.namedtuple('Settings', list(_defaults.keys()))
    
    for k in _defaults.keys():
        
        val = getattr(defaults, k)
        
        setattr(settings, k, val)
    
    globals()['settings'] = settings


def setup(**kwargs):
    
    for param, value in kwargs.items():
        
        setattr(settings, param, value)


def get(param):
    
    if isinstance(param, common.basestring) and hasattr(settings, param):
        
        return getattr(settings, param)


def get_default(param):
    
    if isinstance(param, common.basestring) and hasattr(defaults, param):
        
        return getattr(defaults, param)


def reset(param):
    
    setup(param, get_default(param))


defaults = common._const()

for k, v in _defaults.items():
    
    setattr(defaults, k, v)

reset_all()
