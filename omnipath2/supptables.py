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

import imp

import omnipath2
import omnipath2.table


class SuppTableBase(omnipath2.table.TableBase):
    
    
    def __init__(self, dataset, **kwargs):
        
        self.dataset = dataset
        
        param = {}
        param.update(kwargs)
        
        omnipath2.table.TableBase.__init__(self, **param)
    
    
    def load(self):
        
        self.load_summaries()
        self.load_table()
    
    
    def load_summaries(self):
        
        self.database = omnipath2.data.get_db(self.dataset)
        self.database.update_summaries()
        self.data = self.database.summaries_tab(return_table = True)
    
    
    def load_table(self):
        
        self.header = self.data.pop(0)


class NetworkS2(SuppTableBase):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'fname': 'network_s2_tsv',
            'fname_param': (
                kwargs['dataset'],
                '',
            ),
        }
        param.update(kwargs)
        
        SuppTableBase.__init__(self, **param)
        
        self.export_with_compact_header()
    
    
    def export_with_compact_header(self):
        
        data = self.database.summaries_tab(
            return_table = True,
            label_type = 0,
        )
        header = list(data.pop(0))
        
        tbl = omnipath2.table.TableBase(
            fname = 'network_s2_tsv',
            fname_param = (
                self.dataset,
                '_for-r',
            ),
            data = self.data,
            header = header,
        )


class NetworkS2_PPIall(NetworkS2):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'dataset': 'omnipath',
        }
        param.update(kwargs)
        
        NetworkS2.__init__(self, **param)


class NetworkS2_PPIcurated(NetworkS2):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'dataset': 'curated',
        }
        param.update(kwargs)
        
        NetworkS2.__init__(self, **param)


class NetworkS2_TFtarget(NetworkS2):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'dataset': 'tf_target',
        }
        param.update(kwargs)
        
        NetworkS2.__init__(self, **param)


class NetworkS2_miRNAmRNA(NetworkS2):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'dataset': 'mirna_mrna',
        }
        param.update(kwargs)
        
        NetworkS2.__init__(self, **param)


class NetworkS2_TFmiRNA(NetworkS2):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'dataset': 'tf_mirna',
        }
        param.update(kwargs)
        
        NetworkS2.__init__(self, **param)


class EnzSubS3(SuppTableBase):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'fname': 'enzsub_s3_tsv',
            'dataset': 'enz_sub',
        }
        param.update(kwargs)
        
        SuppTableBase.__init__(self, **param)


class ComplexesS4(SuppTableBase):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'fname': 'complexes_s4_tsv',
            'dataset': 'complex',
        }
        param.update(kwargs)
        
        SuppTableBase.__init__(self, **param)


class AnnotationsS5(SuppTableBase):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'fname': 'annotations_s5_tsv',
            'dataset': 'annotations',
        }
        param.update(kwargs)
        
        SuppTableBase.__init__(self, **param)


class IntercellS6(SuppTableBase):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'fname': 'intercell_s6_tsv',
            'dataset': 'intercell',
        }
        param.update(kwargs)
        
        SuppTableBase.__init__(self, **param)
