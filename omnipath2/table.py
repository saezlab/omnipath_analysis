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

from future.utils import iteritems

import itertools

import pandas as pd

from pypath import session_mod

import omnipath2.path


class TableBase(omnipath2.path.PathBase):
    
    
    def __init__(
            self,
            fname = None,
            fname_param = (),
            fname_timestamp = True,
            timestamp_override = None,
            filetype = 'tsv',
            tables_dir = None,
            dir_timestamp = True,
            data = None,
            sep = '\t',
            header = True,
            log_label = None,
            **kwargs
        ):
        
        if not hasattr(self, '_logger'):
            
            session_mod.Logger.__init__(self, name = log_label or 'op2.table')
        
        tables_dir = tables_dir or omnipath2.data.tables_dir
        
        for k, v in itertools.chain(iteritems(locals()), iteritems(kwargs)):
            
            # we check this because derived classes might have set
            # already attributes
            if not hasattr(self, k) or getattr(self, k) is None:
                
                setattr(self, k, v)
        
        omnipath2.path.PathBase.__init__(
            self,
            fname = self.fname,
            fname_param = self.fname_param,
            fname_timestamp = self.fname_timestamp,
            timestamp_override = self.timestamp_override,
            filetype = self.filetype,
            target_dir = self.tables_dir,
        )
    
    
    def main(self):
        
        omnipath2.path.PathBase.main(self)
        self.load()
        self.export()
    
    
    def load(self):
        
        self.data = self.data or []
        self.header = []
    
    
    def export(self, fname = None):
        
        path = fname or self.path
        
        if isinstance(self.data, pd.DataFrame):
            
            self.data.to_csv(
                path,
                sep = self.sep,
                index = False,
            )
            return
        
        with open(path, 'w') as fp:
            
            if self.header is not None:
                
                _ = fp.write(self.sep.join(self.header))
                _ = fp.write('\n')
            
            fp.write(
                '\n'.join(
                    self.sep.join(
                        str(field) for field in line
                    )
                    for line in self.data
                )
            )
    
    
    def to_data_frame(self):
        """
        Converts the data to ``pandas.DataFrame``.
        """
        
        if not isinstance(self.data, pd.DataFrame):
            
            self.data = pd.DataFrame(
                self.data,
                columns = self.header,
            )
