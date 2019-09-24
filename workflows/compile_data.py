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
    # pickles
   'omnipath_pickle': 'omnipath_pw_es_lr_20190918.pickle',
   'curated_pickle': 'network_curated_20190924.pickle',
   'complex_pickle': 'complexes.pickle',
   'annotations_pickle': 'annotations2.pickle',
   'intercell_pickle': 'intercell.pickle',
   'enz_sub_pickle': 'ptms.pickle',
   'tfregulons_pickle': 'tfregulons.pickle',
   'tfmirna_pickle': 'tfmirna.pickle',
   'mirna_mrna_pickle': 'mirna_mrna.pickle',
   
   # tables
   'network_complete_tsv': 'network_summaries.tsv',
   'network_curated_tsv': 'network-curated_summaries.tsv',
   'network_tf_tsv': 'network-tf.tsv',
   'network_mirna_target_tsv': 'network-mirna-target_summaries.tsv',
   'network_tf_mirna_tsv': 'network-tf-mirna_summaries.tsv',
   'annotations_tsv': 'annotations.tsv',
   'enzsub_tsv': 'enzsub_summaries.tsv',
   'intercell_tsv': 'intercell_summaries.tsv',
   'complex_tsv': 'complex_summaries.tsv',
   
   # tfregulons levels
   'tfregulons_levels': {'A', 'B', 'C', 'D'},
}


class Database(session_mod.Logger)
    
    
    def __init__(self):
        
        session_mod.Logger.__init__(self, name = 'omnipath2.database')
        
        for key, val in itertools.chain(defaults, kwargs.items()):
            
            setattr(self, key, val)
    
    
    def build(self):
        
        self.build_network_complete()
        self.build_network_curated()
        self.build_complex()
    
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
