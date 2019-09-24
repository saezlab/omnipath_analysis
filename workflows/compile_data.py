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

import os
import sys
import imp
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

import workflows.settings as op2_settings


class Database(session_mod.Logger):
    
    
    def __init__(self, rebuild = False, **kwargs):
        
        session_mod.Logger.__init__(self, name = 'op2.database')
        
        self.param = kwargs
        self.rebuild = rebuild
        self.datasets = self.get_param('datasets')
        self.ensure_dirs()
        
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
    
    
    def ensure_dataset(self, dataset):
        
        rebuild_dataset = self.get_param('rebuild_%s' % dataset)
        
        if (
            self.rebuild or
            rebuild_dataset or
            not self.pickle_exists(dataset)
        ):
            
            self.build_dataset(dataset)
            
        else:
            
            self.load_dataset(dataset)
    
    
    def ensure_dirs(self):
        
        if self.get_param('timestamp_dirs'):
            
            self.tables_dir = os.path.join(
                self.get_param('tables_dir'),
                self.timestamp()
            )
            self.figures_dir = os.path.join(
                self.get_param('figures_dir'),
                self.timestamp(),
            )
        
        for _dir in ('pickle', 'tables', 'figures'):
            
            path = self.get_param('%s_dir' % _dir)
            os.makedirs(path, exist_ok = True)
            self._log('%s directory: `%s`.' % (_dir.capitalize(), path))
    
    
    def pickle_path(self, dataset):
        
        return os.path.join(
            self.get_param('pickle_dir'),
            self.get_param('%s_pickle' % dataset),
        )
    
    
    def pickle_exists(self, dataset):
        
        return os.path.exists(self.get_param('%s_pickle' % dataset))
    
    
    def table_path(dataset):
        
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
    
    
    @staticmethod
    def timestamp():
        
        return time.strftime('%Y%m%d')
    
    
    def get_args_curated(self):
        
        resources = copy.deepcopy(data_formats.pathway)
        resources.update(copy.deepcopy(data_formats.enzyme_substrate))
        
        return {'lst': resources}
    
    
    def get_args_tf_target(self):
        
        transcription = copy.deepcopy(data_formats.transcription)
        transcription['tfregulons'].input_args = {
            'levels': self.tfregulons_levels,
        }
        
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
        
        
        delattr(self, dataset)
    
    
    def remove_all(self):
        
        self.foreach_dataset(method = self.ensure_module)
        self.foreach_dataset(method = self.remove_db)
    
    
    def get_param(self, key):
        
        if key in self.param:
            
            return self.param[key]
        
        return op2_settings.get(key)

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
