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
            etattr(self, '%s_pickle' % dataset),
        )
    
    
    def pickle_exists(self, dataset):
        
        return os.path.exists(self.pickle_path(dataset))
    
    
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
        
        if hasattr(self, 'get_%s_args'):
            
            args.update(getattr(self, 'get_%s_args' % dataset)())
        
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
    
    
    def build_network_complete(self):
        
        main.init_db(
            use_omnipath = True,
            kinase_substrate_extra = True,
            ligand_receptor_extra = True,
            pathway_extra = True,
        )
        main.db.save_network(pfile = self.omnipath_pickle)
    
    
    def build_network_curated(self):
        
        curated = main.PyPath()
        resources = copy.deepcopy(data_formats.pathway)
        resources.update(copy.deepcopy(data_formats.enzyme_substrate))
        curated.init_network(resources)
        curated.save_network(pfile = self.curated)
    
    
    def build_complex(self):
        
        co = complex.get_db()
        co.save_to_pickle(self.complex_pickle)
    
    
    def build_annot(self):
        
        a = annot.get_db()
        a.save_to_pickle(annotation_pickle)
    
    
    def build_intercell(self):
        
        i = intercell.IntercellAnnotation()
        i.save_to_pickle(self.intercell_pickle)
    
    
    def build_enzyme_substrate(self):
        
        m = ptm.get_db()
        m.save_to_pickle(self.enz_sub_pickle)
    
    
    def build_transcription(self):
        
        transcription = copy.deepcopy(data_formats.transcription)
        transcription['tfregulons'].input_args = {
            'levels': self.tfregulons_levels,
        }
        
        tr = main.PyPath()
        tr.init_network(transcription)
        tr.save_network(pfile = self.tfregulons_pickle)
    
    
    def build_tf_mirna(self):
        
        tfm = main.PyPath()
        tfm.init_network(data_formats.tf_mirna)
        tfm.save_network(pfile = self.tfmirna_pickle)
    
    
    def build_mirna_target(self):
        
        mit = main.PyPath()
        mit.init_network(data_formats.mirna_target)
        mit.save_network(pfile = self.mirna_mrna_pickle)
    
    
    def load(self):
        
        self.cplex = complex.get_db(pickle_file = self.complex_pickle)
        self.annot = annot.get_db(pickle_file = self.annotations_pickle)
        self.intercell = intercell.get_db(pickle_file = self.intercell_pickle)
        self.omnipath = main.get_db(pfile = self.omnipath_pickle)
        self.curated = main.PyPath()
        self.curated.init_network(pfile = self.curated_pickle)
        self.enz_sub = ptm.get_db(pickle_file = self.enz_sub_pickle)
        self.trascr = main.PyPath()
        self.transcr.init_network(pfile = self.tfregulons_pickle)
        self.tf_mirna = main.PyPath()
        self.tf_mirna.init_network(pfile = self.tfmirna_pickle)
        self.mirna_mrna = main.PyPath()
        self.mirna_mrna.init_network(pfile = self.mirna_mrna_pickle)
    
    
    def compile_tables(self):
        
        
        self.omnipath.update_summaries()
        self.omnipath.summaries_tab(outfile = self.network_complete_tsv)
        
        self.omnipath.update_summaries()
        self.omnipath.summaries_tab(outfile = self.network_curated_tsv)
        
        self.transcr.update_summaries()
        self.transcr.summaries_tab(outfile = self.network_tf_tsv)
        
        self.tf_mirna.update_summaries()
        self.tf_mirna.summaries_tab(outfile = self.network_tf_mirna_tsv)
        
        self.mirna_mrna.update_summaries()
        self.mirna_mrna.summaries_tab(outfile = self.network_mirna_target_tsv)
        
        self.cplex.update_summaries()
        self.cplex.summaries_tab(self.complex_tsv)
        
        self.annot.update_summaries()
        self.annot.summaries_tab(outfile = self.annotations_tsv)
        
        self.intercell.update_summaries()
        self.summaries_tab(outfile = self.intercell_tsv)

# count annotation records

# proteins
sum([aaa.numof_protein_records() for aaa in a.annots.values()]) # 1972327

# complexes
sum([aaa.numof_complex_records() for aaa in a.annots.values()]) # 218360

# miRNAs
sum([aaa.numof_mirna_records() for aaa in a.annots.values()]) # 355

# total
sum([aaa.numof_records() for aaa in a.annots.values()]) # 2191042

# count annotated entities

all_annotated_proteins = set.union(
    *(set(aaa.all_proteins()) for aaa in a.annots.values())
)

all_uniprots = set(uniprot_input.all_uniprots(swissprot = 'yes'))

len(all_annotated_proteins & all_uniprots) # 20066

# number of annotated complexes
all_annotated_complexes = set.union(
    *(set(aaa.all_complexes()) for aaa in a.annots.values())
)

len(all_annotated_complexes) # 16477

# number of complexes in OmniPath:
all_complexes = set(complex.db.complexes.values())
len(all_complexes)


cplex_stats = c.summaries_tab(outfile = 'complex_stats.tsv')
