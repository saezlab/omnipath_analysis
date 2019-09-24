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

import imp
import pprint
import copy
import collections
import itertools

from pypath import mapping
from pypath import dataio
from pypath import data_formats
from pypath import intercell_annot
from pypath import annot
from pypath import intercell
from pypath import main
from pypath import network
from pypath import uniprot_input
from pypath import complex
from pypath import ptm
from pypath import session_mod

import settings as op2_settings


class Database(session_mod.Logger)
    
    
    def __init__(self, rebuild = False):
        
        session_mod.Logger.__init__(self, name = 'op2.database')
        
        self.param = kwargs
        self.rebuild = rebuild
        self.datasets = op2_settings.get('datasets')
        self.ensure_dirs()
        
        self._log('OmniPath2 database builder initialized.')
    
    
    def build(self):
        
        self._log(
            'Building databases. Rebuild forced: %s.' % str(self.rebuild)
        )
        
        self.foreach_dataset(method = self.ensure_dataset)
    
    
    def ensure_dataset(self, dataset):
        
        rebuild_dataset = op2_settings.get('rebuild_%s' % dataset)
        
        if (
            self.rebuild or
            rebuild_dataset or
            not self.pickle_exists(dataset)
        ):
            
            self.build_dataset(dataset)
            
        else:
            
            self.load_dataset(dataset)
    
    
    def ensure_dirs(self):
        
        if op2_settings.get('timestamp_dirs'):
            
            self.tables_dir = os.path.join(
                op2_settings.get('tables_dir'),
                self.timestamp()
            )
            self.figures_dir = os.path.join(
                op2_settings.get('figures_dir'),
                self.timestamp(),
            )
        
        for _dir in ('pickle', 'tables', 'figures'):
            
            path = op2_settings.get('%s_dir' % _dir)
            os.makedirs(path, exist_ok = True)
            self._log('%s directory: `%s`.' % (_dir.capitalize(), path))
    
    
    def pickle_path(self, dataset):
        
        return os.path.join(
            op2_settings.get('pickle_dir'),
            getattr(self, '%s_pickle' % dataset),
        )
    
    
    def pickle_exists(self, dataset):
        
        return os.path.exists(op2_settings.get('dataset'))
    
    
    def table_path(dataset):
        
        return os.path.join(
            op2_settings.get('tables_dir'),
            op2_settings.get('%s_tsv' % dataset),
        )
    
    
    def build_dataset(self, dataset):
        
        self._log('Building dataset `%s`.' % dataset)
        
        args = get_build_args(dataset)
        
        mod = self.ensure_module(dataset)
        
        db = mod.get_db(**args)
        
        pickle_path = self.pickle_path(dataset)
        self._log('Saving dataset `%s` to `%s`.' % (dataset, pickle_path))
        db.save_to_pickle(pickle_file = pickle_path)
        
        self._log('Successfully built dataset `%s`.' % dataset)
        
        setattr(self, dataset, db)
    
    
    def ensure_module(self, dataset):
        
        mod = op2_settings.get('%s_mod' % dataset)
        
        if hasattr(mod, 'db'):
            
            delattr(mod, 'db')
        
        return mod
    
    
    def get_build_args(self, dataset):
        
        args = op2_settings.get('%s_args' % dataset) or {}
        
        if hasattr(self, 'get_args_%s' % dataset):
            
            args.update(getattr(self, 'get_args_%s' % dataset)())
        
        return args
    
    
    def load_dataset(self, dataset):
        
        pickle_path = self.pickle_path(dataset)
        
        self._log('Loading dataset `%s` from `%s`.' % (dataset, pickle_path))
        
        mod = self.ensure_module(dataset)
        
        mod.get_db(pickle_file = pickle_path)
        
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
            
            method(dataset)
    
    
    def get_db(self, dataset):
        
        self.ensure_dataset(dataset)
        
        return getattr(self, dataset)
    
    
    def remove_db(self, dataset):
        
        delattr(self, dataset)

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
