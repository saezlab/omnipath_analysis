#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Denes Turei 2019
# turei.denes@gmail.com

from workflows import figures_preprocess as figpreproc

omnipath_pickle = 'omnipath.pickle'
annotation_pickle = 'annotations2.pickle'
complex_pickle = 'complexes.pickle'
intercell_pickle = 'intercell.pickle'

proc = figpreproc.FiguresPreprocess(
    omnipath_pickle = omnipath_pickle,
    annotation_pickle = annotation_pickle,
    complex_pickle = complex_pickle,
)

proc.setup()
proc.load()

