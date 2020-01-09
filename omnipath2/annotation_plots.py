#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019-2020 Saez Lab
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

import itertools

import pattern.en

import pypath.common as common
import pypath.session_mod as session_mod
import pypath.annot as annot

import omnipath2
import omnipath2.settings as op2_settings
import omnipath2.plot as plot
import omnipath2.intercell_plots


class EntitiesByResource(omnipath2.intercell_plots.CountsScatterBase):
    
    
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
        
        omnipath2.intercell_plots.CountsScatterBase.__init__(self, **param)
    
    
    def load_data(self):
        
        self.annot = omnipath2.data.get_db('annotations')
        self.counts = dict(
            (
                resource,
                annot._numof_entities(entity_types = self.entity_types)
            )
            for resource, annot in self.annot.annots.items()
        )
        
        omnipath2.intercell_plots.CountsScatterBase.load_data(self)


class RecordsByResource(omnipath2.intercell_plots.CountsScatterBase):
    
    
    def __init__(
            self,
            entity_types = 'protein',
            **kwargs
        ):
        
        self.entity_types = common.to_set(entity_types)
        
        param = {
            'fname': 'annot_records_by_resource_pdf',
            'fname_param': (
                '-'.join(et for et in self.entity_types)
                    if self.entity_types else
                'all'
            ),
            'title': 'Number of %s records\nby annotation resource' % (
                entity_types
                    if isinstance(entity_types, common.basestring) else
                'all'
            ),
            'xlab': 'Number of %s records' % (
                entity_types
                    if isinstance(entity_types, common.basestring) else
                'all'
            ),
            'ylab': 'Annotation resources',
            'height': 7,
            'xscale_log': True,
            'xlim': (1, 1e7),
        }
        param.update(kwargs)
        
        omnipath2.intercell_plots.CountsScatterBase.__init__(self, **param)
    
    
    def load_data(self):
        
        self.annot = omnipath2.data.get_db('annotations')
        self.counts = dict(
            (
                resource,
                _annot.numof_records(entity_types = self.entity_types)
            )
            for resource, _annot in self.annot.annots.items()
        )
        
        omnipath2.intercell_plots.CountsScatterBase.load_data(self)


class AnnotationNetworkOverlap(omnipath2.intercell_plots.CountsScatterBase):
    
    
    def __init__(
            self,
            entity_types = 'protein',
            network_dataset = 'omnipath',
            **kwargs
        ):
        
        self.entity_types = common.to_set(entity_types)
        self.network_dataset = network_dataset
        
        param = {
            'fname': 'annot_entities_in_network_pdf',
            'fname_param': (
                '-'.join(et for et in self.entity_types)
                    if self.entity_types else
                'all',
                self.network_dataset,
            ),
            'title': ' %s in the network (%s)\nby annotation resource' % (
                pattern.en.pluralize(entity_types).capitalize()
                    if isinstance(entity_types, common.basestring) else
                'all',
                self.network_dataset,
            ),
            'xlab': ' %s in the network' % (
                pattern.en.pluralize(entity_types).capitalize()
                    if isinstance(entity_types, common.basestring) else
                'Entities'
            ),
            'ylab': 'Annotation resources',
            'height': 7,
            'xscale_log': True,
            'xlim': (1, 1e5),
        }
        param.update(kwargs)
        
        omnipath2.intercell_plots.CountsScatterBase.__init__(self, **param)
    
    
    def load_data(self):
        
        self.network = omnipath2.data.get_db(self.network_dataset)
        self.annot = omnipath2.data.get_db('annotations')
        self.in_network = {
            _id
            for _id in self.network.get_identifiers(
                entity_type = self.entity_types
            )
        }
        self.counts = dict(
            (
                resource,
                len(
                    set(annot.all_entities(self.entity_types)) &
                    self.in_network
                )
            )
            for resource, annot in self.annot.annots.items()
        )
        
        omnipath2.intercell_plots.CountsScatterBase.load_data(self)
