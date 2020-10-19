#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019-2020 Saez Lab
#
# OmniPath2 analysis and figures suite
#
# Authors:
#
# Dénes Türei
# turei.denes@gmail.com
#
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://omnipathdb.org/
#

import os

import pypath.omnipath as omnipath
import pypath.core.network as network
import pypath.resources.network as netres
import pypath.omnipath.export as export


def export_network(inputdef, name):

    netw = (
        inputdef
            if isinstance(inputdef, network.Network) else
        network.Network(inputdef)
    )

    outfile = os.path.join('data', '%s.tsv' % name)
    exp = export.Export(
        network = netw,
        outfile = outfile,
    )
    exp.webservice_interactions_df()
    exp.write_tab()




export_network(netres.pathwaycommons, 'ppi_pc')
export_network(netres.pathwaycommons_transcription, 'tftarget_pc')
export_network(omnipath.db.get_db('curated'), 'ppi_omnipath-curated')

