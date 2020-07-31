#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019 Saez Lab
#
# OmniPath2 analysis and figures suite
#
# Authors:
#
# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de
#
# Dénes Türei
# turei.denes@gmail.com
#
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://omnipathdb.org/
#

import importlib as imp

import omnipath2.main
from omnipath2 import settings
from pypath.resources import data_formats
from pypath.resources import network as netres
from pypath.core import annot

def kegg_off():

    settings.setup(tf_target_pickle = 'tftarget_paper.pickle')
    settings.setup(omnipath_pickle = 'omnipath_paper.pickle')
    settings.setup(annotations_pickle = 'annotations_paper.pickle')
    del data_formats.pathway_noref['kegg']
    del data_formats.transcription_onebyone['kegg']
    del data_formats.transcription['kegg']
    annot.protein_sources_default.discard('KEGG')
    imp.reload(netres)



kegg_off()

workflow = omnipath2.main.Main()

workflow.main()
