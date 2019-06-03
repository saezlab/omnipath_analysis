#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

import imp
import pprint
import collections

from pypath import mapping
from pypath import dataio
from pypath import intercell_annot
from pypath import annot
from pypath import intercell
from pypath import main
from pypath import network


omnipath_pickle = 'omnipath_no-complex-expansion_20190603.pickle'
complex_pickle = 'complexes.pickle'
annotations_pickle = 'annotations.pickle'
intercell_pickle = 'intercell.pickle'


pa = main.PyPath()
pa.load_omnipath(
    kinase_substrate_extra = True,
    ligand_receptor_extra = True,
)

pa.save_network(pfile = omnipath_pickle)


co = complexes.get_db()
co.save_to_pickle(complex_pickle)


a = annot.get_db()
a.save_to_pickle(annotation_pickle)


i = intercell.IntercellAnnotation()
i.save_to_pickle(intercell_pickle)
