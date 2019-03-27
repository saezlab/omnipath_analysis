#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

import imp

from pypath import mapping
from pypath import dataio
from pypath import complex


def reload():
    
    imp.reload(mapping)
    mapping.init()
    imp.reload(dataio)
    imp.reload(complex)


co = complex.ComplexAggregator()

# this is a dict of Complex objects:
co.complexes

# with stoichiometries:
co.complexes['O15151-P0CG47-P62837-Q00987'].components

# sources and references:
co.complexes['P01584-P27930-Q9NPH3'].sources
co.complexes['P18507-P28472-Q16445'].references
