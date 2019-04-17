#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

import imp
import pprint

from pypath import mapping
from pypath import dataio
from pypath import intercell_annot
from pypath import annot
from pypath import intercell


def reload():
    
    imp.reload(intercell_annot)
    imp.reload(annot)
    imp.reload(intercell)


# this attempts to come to a consensus about the roles of the proteins
# then the next class will check these against the network resources
# and provide reports for plotting
i = intercell.IntercellAnnotation()

pprint.pprint([(k, len(v)) for k, v in i.classes.items()])

i.export(
    fname = 'intercell_classes_20190411.tsv',
    sep = '\t',
    index = False,
)
