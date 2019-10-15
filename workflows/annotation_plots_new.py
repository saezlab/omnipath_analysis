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

import itertools

import pattern.en

import pypath.common as common

import workflows
import workflows.settings as op2_settings
import workflows.plot as plot
import workflows.intercell_plots_new as intercell_plots


class EntitiesByResource(intercell_plots.CountsScatterBase):
    
    
    def __init__(
            self,
            entity_types = 'protein',
            **kwargs
        ):
        
        self.entity_types = common.to_set(entity_types)
        
        param = {
            'fname': 'annot_entities_by_resource_pdf',
            'fname_param': (
                '-'.join(pattern.en.pluralize(et) for et in self.entity_types)
                    if self.entity_types else
                'all-entities'
            ),
            'title': 'Number of %s by annotation resource' % (
                pattern.en.pluralize(entity_types)
                    if isinstance(entity_types, common.basestring) else
                'entities'
            ),
            'xlab': 'Number of %s' % (
                pattern.en.pluralize(entity_types)
                    if isinstance(entity_types, common.basestring) else
                'entities'
            ),
            'ylab': 'Annotation resources',
            'height': 7,
            'xscale_log': True,
            'xlim': (1, 1e5),
        }
        param.update(kwargs)
        
        intercell_plots.CountsScatterBase.__init__(self, **param)
    
    
    def load_data(self):
        
        self.annot = workflows.data.get_db('annotations')
        self.counts = dict(
            (
                resource,
                annot._numof_entities(entity_types = self.entity_types)
            )
            for resource, annot in self.annot.annots.items()
        )
        
        intercell_plots.CountsScatterBase.load_data(self)
