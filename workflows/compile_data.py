#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

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


defaults = {
    # pickle_dir
    'pickle_dir': 'pickles',
    
    # tables dir
    'tables_dir': 'tables',
    
    # pickles
    'omnipath_pickle': 'omnipath_pw_es_lr_20190918.pickle',
    'curated_pickle': 'network_curated_20190924.pickle',
    'complex_pickle': 'complexes.pickle',
    'annotations_pickle': 'annotations2.pickle',
    'intercell_pickle': 'intercell.pickle',
    'enz_sub_pickle': 'ptms.pickle',
    'tf_target_pickle': 'tfregulons.pickle',
    'tf_mirna_pickle': 'tfmirna.pickle',
    'mirna_mrna_pickle': 'mirna_mrna.pickle',

    # tables
    'omnipath_tsv': 'network_summaries.tsv',
    'curated_tsv': 'network-curated_summaries.tsv',
    'tf_target_tsv': 'network-tf.tsv',
    'mirna_mrna_tsv': 'network-mirna-target_summaries.tsv',
    'tf_mirna_tsv': 'network-tf-mirna_summaries.tsv',
    'annotations_tsv': 'annotations.tsv',
    'enz_sub_tsv': 'enzsub_summaries.tsv',
    'intercell_tsv': 'intercell_summaries.tsv',
    'complex_tsv': 'complex_summaries.tsv',

    # tfregulons levels
    'tfregulons_levels': {'A', 'B', 'C', 'D'},

    # datasets
    [
       'omnipath',
       'curated',
       'complex',
       'annotations',
       'intercell',
       'tf_target',
       'tf_mirna',
       'mirna_target',
       'enz_sub',
    ],
    
    'timestamp_dirs': True,
    
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
}


class Database(session_mod.Logger)
    
    
    def __init__(self, rebuild = False):
        
        session_mod.Logger.__init__(self, name = 'omnipath2.database')
        
        for key, val in itertools.chain(defaults, kwargs.items()):
            
            setattr(self, key, val)
        
        self.rebuild = rebuild
        
        self._log('OmniPath2 database builder initialized.')
    
    
    def main(self):
        
        self.ensure_dirs()
        self.build()
        build()
        load()
    
    
    def build(self):
        
        self._log(
            'Building databases. Rebuild forced: %s.' % str(self.rebuild)
        )
        
        for dataset in self.datasets:
            
            self.ensure_dataset(dataset)
    
    
    def ensure_dataset(self, dataset):
        
        rebuild_attr = 'rebuild_%s' % dataset
        
        if (
            self.rebuild or
            not self.pickle_exists(dataset) or (
                hasattr(self, rebuild_attr) and
                getattr(self, rebuild_attr)
            )
        ):
            
            self.build_dataset(dataset)
            
        else:
            
            self.load_dataset(dataset)
    
    
    def ensure_dirs(self):
        
        if self.timestamp_dirs:
            
            self.tables_dir = os.path.join(
                self.tables_dir,
                self.timestamp()
            )
            self.figures_dir = os.path.join(
                self.figures_dir,
                self.timestamp(),
            )
        
        for _dir in ('pickle', 'tables', 'figures'):
            
            path = getattr(self, '%s_dir' % _dir
            os.makedirs(path, exist_ok = True)
            self._log('%s directory: `%s`.' % (_dir.capitalize(), path))
    
    
    def pickle_path(self, dataset):
        
        return os.path.join(
            self.pickle_dir,
            getattr(self, '%s_pickle' % dataset),
        )
    
    
    def pickle_exists(self, dataset):
        
        return os.path.exists(self.pickle_path(dataset))
    
    
    def table_path(self, dataset):
        
        return os.path.join(
            self.tables_dir,
            getattr(self, '%s_tsv' % dataset),
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
        
        mod = getattr(self, '%s_mod' % dataset)
        
        if hasattr(mod, 'db'):
            
            delattr(mod, 'db')
        
        return mod
    
    
    def get_build_args(self, dataset):
        
        args = (
            getattr(self, '%s_args' % dataset)
                if hasattr(self, '%s_args' % dataset) else
            {}
        )
        
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
        
        for dataset in self.datasets:
            
            self.compile_table(dataset)
    
    
    def compile_table(self, dataset):
        
        table_path = self.table_path(dataset)
        db = self.get_db(dataset)
        db.update_summaries()
        db.summaries_tab(outfile = table_path)
    
    
    def get_db(self, dataset):
        
        return getattr(self, dataset)

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
