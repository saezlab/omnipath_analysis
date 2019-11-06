#!/usr/bin/env python
#-*- coding: utf-8 -*-

#
#  This file is part of the `omnipath2` Python module
#
#  Copyright
#  2019
#  Heidelberg University, Uniklinik RWTH Aachen
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections
import itertools

import omnipath2
import omnipath2.plot
import omnipath2.intercell_plots_new as intercell_plots


class ComplexesByResource(intercell_plots.CountsScatterBase):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'fname': 'complexes_by_resource_pdf',
            'title': 'Number of complexes by resource',
            'xlab': 'Number of complexes',
            'ylab': 'Complex resources',
            'height': 4,
            'xscale_log': True,
            'xlim': (1, 1e5),
        }
        param.update(kwargs)
        
        intercell_plots.CountsScatterBase.__init__(self, **param)
    
    
    def load_data(self):
        
        self.complexes = omnipath2.data.get_db('complex')
        self.counts = collections.Counter(
            itertools.chain(*(
                cplex.sources
                for cplex in self.complexes.complexes.values()
            ))
        )
        
        intercell_plots.CountsScatterBase.load_data(self)


