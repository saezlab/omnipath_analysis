#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019-2020 Saez Lab
#
# OmniPath2 analysis and figures suite
#
# Authors:
#
# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de
#
# Dénes Türei
# turei.denes@gmail.com
#
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://omnipathdb.org/
#

import collections

from pypath.share import common


_defaults = {

    'files_json': 'files.json',

    # pickle dumps of all databases
    'pickle_dir': 'pickles',

    # directory for exported tables
    'tables_dir': 'tables',

    # directory for figures
    'figures_dir': 'figures',

    # directory for LaTeX
    'latex_dir': 'latex',

    # directory for auto-generated LaTeX snippets
    'latex_auto_dir': 'inc_auto',

    # directory with palettes
    'palette_dir': 'palettes',

    'timestamp_dirs': True,

    'timestamp_format': '%Y%m%d',

    # pickles
    'omnipath_pickle': 'network_omnipath.pickle',
    'curated_pickle': 'network_curated.pickle',
    'complex_pickle': 'complexes.pickle',
    'annotations_pickle': 'annotations.pickle',
    'intercell_pickle': 'intercell.pickle',
    'enz_sub_pickle': 'enz_sub_9606.pickle',
    'tf_target_pickle': 'tftarget.pickle',
    'tf_mirna_pickle': 'tfmirna.pickle',
    'mirna_mrna_pickle': 'mirna_mrna.pickle',
    'lncrna_mrna_pickle': 'lncrna_mrna.pickle',

    # supplementary tables
    'network_s2_tsv': 'S2_network_%s%s',
    'enzsub_s3_tsv': 'S3_enz-sub',
    'complexes_s4_tsv': 'S4_complexes',
    'annotations_s5_tsv': 'S5_annotations',
    'intercell_s6_tsv': 'S6_intercell',

    # for R plotting
    'intercell_classes_tsv': 'intercell_classes', #
    'main_coverage_tsv': 'main_coverage_%s', #
    'intercell_network_by_resource_tsv': 'stats_by_resource_%s_%s', #
    'intercell_network_tsv': 'intercell_network_by_resource', #
    'resources_by_entity_tsv': 'resources_by_entity_%s', #
    'intercell_annots_by_entity_tsv': 'intercell_annotations_by_entity', #
    'annots_by_entity_tsv': 'annotations_by_entity', #
    'enz_sub_tsv': 'enzyme_substrate', #
    'complexes_by_resource_tsv': 'complexes_by_resource_tsv', #
    'connections_tsv': 'connections_%s_%s', #
    'category_overlaps_tsv': 'category_overlaps', #
    'network_coverage_tsv': 'network_coverage_%s', #
    'inter_class_summary_tsv': 'inter_class_network_summary_%s_%s_%s',
    'network_consistency_tsv': 'network_consistency_%s',

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
       'lncrna_mrna',
       'enz_sub',
    ],

    'omnipath_mod': 'network',
    'curated_mod': 'network',
    'complex_mod': 'complex',
    'annotations_mod': 'annot',
    'intercell_mod': 'intercell',
    'enz_sub_mod': 'enz_sub',
    'tf_target_mod': 'network',
    'tf_mirna_mod': 'network',
    'mirna_mrna_mod': 'network',
    'lncrna_mrna_mod': 'network',

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
    'counts_by_class_pdf': 'counts_by_class_%s',
    'counts_by_resource_pdf': 'counts_by_resource_%s',
    'inter_class_sim_pdf': 'inter_class_sim',
    'inter_class_chordplot_pdf': 'inter_class_chordplot_%s_%s_%s',
    'annot_entities_by_resource_pdf': 'annot_%s_by_resource',
    'annot_records_by_resource_pdf': 'annot_%s_records_by_resource',
    'annot_entities_in_network_pdf': 'annot_%s_in_network_%s',
    'netw_node_edge_counts_pdf': 'nodes_edges_%s_%s',
    'complexes_by_resource_pdf': 'complexes_by_resource',

    # default palette
    'palette': 'Nico_3',

    'r_workflow_command': 'Rscript omnipath2_workflow.r',

    'intercell_main_classes': [
        'adhesion',
        'secreted_enzyme',
        'ecm',
        'ligand',
        'receptor',
        'cell_surface_enzyme',
        'cell_surface_ligand',
        'transporter'
    ],

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


class Settings(object):


    def __init__(self, **kwargs):

        self.__dict__.update(kwargs)


Defaults = collections.namedtuple(
    'Defaults',
    sorted(_defaults.keys()),
)


def reset_all():

    settings = Settings()

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


defaults = Defaults(**_defaults)

reset_all()
