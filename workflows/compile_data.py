#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

import imp
import pprint
import copy
import collections

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


omnipath_pickle = 'omnipath_pw_es_lr_20190918.pickle'
complex_pickle = 'complexes.pickle'
annotations_pickle = 'annotations2.pickle'
intercell_pickle = 'intercell.pickle'
ptms_pickle = 'ptms.pickle'
tfregulons_pickle = 'tfregulons.pickle'
tfmirna_pickle = 'tfmirna.pickle'
mirna_mrna_pickle = 'mirna_mrna.pickle'


def build():
    
    main.init_db(
        use_omnipath = True,
        kinase_substrate_extra = True,
        ligand_receptor_extra = True,
        pathway_extra = True,
    )
    main.db.save_network(pfile = omnipath_pickle)
    
    co = complex.get_db()
    co.save_to_pickle(complex_pickle)
    
    a = annot.get_db()
    a.save_to_pickle(annotation_pickle)
    
    i = intercell.IntercellAnnotation()
    i.save_to_pickle(intercell_pickle)
    
    m = ptm.get_db()
    m.save_to_pickle(ptms_pickle)
    
    transcription = copy.deepcopy(data_formats.transcription)
    transcription['tfregulons'].input_args = {
        'levels': {'A', 'B', 'C', 'D'},
    }
    
    tr = main.PyPath()
    tr.init_network(transcription)
    tr.save_network(pfile = tfregulons_pickle)
    
    tfm = main.PyPath()
    tfm.init_network(data_formats.tf_mirna)
    tfm.save_network(pfile = tfmirna_pickle)
    
    
    mit = main.PyPath()
    mit.init_network(data_formats.mirna_target)
    mit.save_network(pfile = mirna_mrna_pickle)



def get_databases():
    
    c = complex.get_db(pickle_file = complex_pickle)
    a = annot.get_db(pickle_file = annotations_pickle)
    i = intercell.get_db(pickle_file = intercell_pickle)
    p = main.get_db(pfile = omnipath_pickle)
    m = ptm.get_db(pickle_file = ptms_pickle)
    tr = main.PyPath()
    tr.init_network(pfile = tfregulons_pickle)
    tfm = main.PyPath()
    tfm.init_network(pfile = tfmirna_pickle)
    mit = main.PyPath()
    mit.init_network(pfile = mirna_mrna_pickle)
    
    return c, a, i, p, m, tr, tfm, mit


c, a, i, p, m, tr, tfm, mit = get_databases()

p.update_summaries()
_ = p.summaries_tab(outfile = 'ppi_network_summaries.tsv')
tr.update_summaries()
_ = tr.summaries_tab(outfile = 'tf_network_summaries.tsv')
tfm.update_summaries()
_ = tfm.summaries_tab(outfile = 'tf-mirna_network_summaries.tsv')
mit.update_summaries()
_ = mit.summaries_tab(outfile = 'mirna-target_network_summaries.tsv')

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
