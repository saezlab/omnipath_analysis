#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

import imp

from pypath import mapping
from pypath import dataio
from pypath import annot


def reload():

    imp.reload(mapping)
    mapping.init()
    imp.reload(dataio)
    imp.reload(annot)


a = annot.AnnotationTable(keep_annotators = True)
a.load()
# there is a boolean array representing the membership of proteins
# in each category; rows are protein, columns are categories
a.data # the array
a.names # the categories (columns)
a.uniprots # the proteins (rows)

# the individual resources are here:
a.annots

# for example:
m = a.annots['Matrisome']

# the data is a dict with a set of annotations for each protein
# annotations are named tuples:
m.annot
# these objects have methods to query them various ways

# e.g. 'Q9Y264' has
# MatrisomeAnnotation(
#   mainclass='Matrisome-associated',
#   subclass='Secreted Factors',
#   subsubclass=None)
# hence
a.cols[('Matrisome', 'Matrisome-associated', 'Secreted Factors')] # 154
a.rows['Q9Y264'] # 19781
a.data[19781, 154] # True
