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


a = annot.AnnotationTable()
a.load()
