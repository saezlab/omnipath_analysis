#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019 Saez Lab
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

import os
import sys
import importlib as imp
import time
import pprint
import copy
import collections
import itertools

from pypath import data_formats
from pypath import annot
from pypath import intercell
from pypath import main
from pypath import complex
from pypath import ptm
from pypath import session_mod
from pypath import network

import omnipath2.settings as op2_settings


class Database(session_mod.Logger):


    def __init__(self, rebuild = False, **kwargs):

        session_mod.Logger.__init__(self, name = 'op2.database')

        self.timestamp = time.strftime(op2_settings.get('timestamp_format'))
        self.param = kwargs
        self.rebuild = rebuild
        self.datasets = self.get_param('datasets')
        self.ensure_dirs()
        self.network_dfs = {}

        self._log('OmniPath2 database builder initialized.')


    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
        self.foreach_dataset(method = self.reload_module)


    def build(self):

        self._log(
            'Building databases. Rebuild forced: %s.' % str(self.rebuild)
        )

        self.foreach_dataset(method = self.ensure_dataset)


    def ensure_dataset(
            self,
            dataset,
            force_reload = False,
            force_rebuild = False,
        ):

        for dep_dataset in self.dataset_dependencies(dataset):

            self.ensure_dataset(dep_dataset)

        rebuild_dataset = self.get_param('rebuild_%s' % dataset)

        if force_reload or force_rebuild or not hasattr(self, dataset):

            if (
                force_rebuild or
                self.rebuild or
                rebuild_dataset or
                not self.pickle_exists(dataset)
            ):

                self.remove_db(dataset)
                self.build_dataset(dataset)

            elif not hasattr(self, dataset) or force_reload:

                self.load_dataset(dataset)


    def dataset_dependencies(self, dataset):

        deps = self.get_param('dependencies')

        return deps[dataset] if dataset in deps else ()


    def ensure_dirs(self):

        if self.get_param('timestamp_dirs'):

            self.tables_dir = os.path.join(
                self.get_param('tables_dir'),
                self.timestamp
            )
            self.figures_dir = os.path.join(
                self.get_param('figures_dir'),
                self.timestamp,
            )
            op2_settings.setup(
                tables_dir = self.tables_dir,
                figures_dir = self.figures_dir,
            )

        for _dir in ('pickle', 'tables', 'figures'):

            path = self.get_param('%s_dir' % _dir)
            os.makedirs(path, exist_ok = True)
            self._log('%s directory: `%s`.' % (_dir.capitalize(), path))


    def pickle_path(self, dataset):

        pickle_fname = self.get_param('%s_pickle' % dataset)

        return os.path.join(
            self.get_param('pickle_dir'),
            pickle_fname,
        )


    def pickle_exists(self, dataset):

        return os.path.exists(self.pickle_path(dataset))


    def table_path(self, dataset):

        return os.path.join(
            self.get_param('tables_dir'),
            self.get_param('%s_tsv' % dataset),
        )


    def build_dataset(self, dataset):

        self._log('Building dataset `%s`.' % dataset)

        args = self.get_build_args(dataset)

        mod = self.ensure_module(dataset)

        db = mod.get_db(**args)

        pickle_path = self.pickle_path(dataset)
        self._log('Saving dataset `%s` to `%s`.' % (dataset, pickle_path))
        db.save_to_pickle(pickle_file = pickle_path)

        self._log('Successfully built dataset `%s`.' % dataset)

        setattr(self, dataset, db)

        self._add_network_df(dataset)


    def ensure_module(self, dataset, reset = True):

        mod_str = self.get_param('%s_mod' % dataset)
        mod = sys.modules['pypath.%s' % mod_str]

        if reset and hasattr(mod, 'db'):

            delattr(mod, 'db')

        return mod


    def reload_module(self, dataset):

        mod = self.ensure_module(dataset, reset = False)
        imp.reload(mod)

        if hasattr(mod, 'db'):

            mod.db.reload()


    def get_build_args(self, dataset):

        args = self.get_param('%s_args' % dataset) or {}

        if hasattr(self, 'get_args_%s' % dataset):

            args.update(getattr(self, 'get_args_%s' % dataset)())

        return args


    def load_dataset(self, dataset):

        pickle_path = self.pickle_path(dataset)

        self._log('Loading dataset `%s` from `%s`.' % (dataset, pickle_path))

        mod = self.ensure_module(dataset)

        setattr(self, dataset, mod.get_db(pickle_file = pickle_path))

        self._log('Loaded dataset `%s` from `%s`.' % (dataset, pickle_path))

        self._add_network_df(dataset)


    def get_args_curated(self):

        resources = copy.deepcopy(data_formats.pathway)
        resources.update(copy.deepcopy(data_formats.enzyme_substrate))

        return {'lst': resources}


    def get_args_tf_target(self):

        transcription = copy.deepcopy(data_formats.transcription)
        dorothea = {}
        
        for level in self.get_param('tfregulons_levels'):
            
            dorothea['dorothea_%s' % level] = copy.deepcopy(
                transcription['dorothea']
            )
            dorothea['dorothea_%s' % level].name = 'DoRothEA_%s' % level
            dorothea['dorothea_%s' % level].input_args = {'levels': {level}}
        
        del transcription['dorothea']
        transcription.update(dorothea)

        return {'lst': transcription}


    def get_args_tf_mirna(self):

        return {'lst': data_formats.tf_mirna}


    def get_args_mirna_mrna(self):

        return {'lst': data_formats.mirna_target}


    def compile_tables(self):

        self.foreach_dataset(method = self.compile_table)


    def compile_table(self, dataset):

        table_path = self.table_path(dataset)
        db = self.get_db(dataset)
        db.update_summaries()
        db.summaries_tab(outfile = table_path)


    def foreach_dataset(self, method):

        for dataset in self.datasets:

            _ = method(dataset)


    def get_db(self, dataset):

        self.ensure_dataset(dataset)

        return getattr(self, dataset)


    def remove_db(self, dataset):

        if hasattr(self, dataset):

            delattr(self, dataset)


    def remove_all(self):

        self.foreach_dataset(method = self.ensure_module)
        self.foreach_dataset(method = self.remove_db)


    def get_param(self, key):

        if key in self.param:

            return self.param[key]

        return op2_settings.get(key)


    def _create_network_df(self, dataset = 'omnipath', **kwargs):

        graph = self.get_db(dataset)

        return self._network_df(graph, **kwargs)

    def network_df(self, dataset, by_source = False):

        self.ensure_dataset(dataset)

        by_source_str = 'by_source' if by_source else 'plain'

        return self.network_dfs[dataset][by_source_str]


    def network_df_by_source(self, dataset = 'omnipath'):

        self.ensure_dataset(dataset)

        return self.network_dfs[dataset]['by_source']


    def _network_df(self, pp_obj, **kwargs):

        return network.Network.from_igraph(pp_obj, **kwargs)


    def _add_network_df(self, dataset):

        obj = getattr(self, dataset)

        if hasattr(obj, 'graph') and isinstance(obj.graph, main.igraph.Graph):

            network_df = self._network_df(obj, by_source = False)
            network_df_by_source = self._network_df(obj, by_source = True)

            self.network_dfs[dataset] = {}
            self.network_dfs[dataset]['plain'] = network_df
            self.network_dfs[dataset]['by_source'] = network_df_by_source

            self._log('Created network data frames for `%s`.' % dataset)


    def set_network(self, dataset, by_source = False, **kwargs):
        """
        Sets dataset as the default
        """

        network_df = self.network_df(dataset, by_source = by_source, **kwargs)

        self.ensure_dataset('intercell')

        self.intercell.register_network(network_df)

#
# to be removed once we have it elsewhere:
#

## count annotation records

## proteins
#sum([aaa.numof_protein_records() for aaa in a.annots.values()]) # 1972327

## complexes
#sum([aaa.numof_complex_records() for aaa in a.annots.values()]) # 218360

## miRNAs
#sum([aaa.numof_mirna_records() for aaa in a.annots.values()]) # 355

## total
#sum([aaa.numof_records() for aaa in a.annots.values()]) # 2191042
