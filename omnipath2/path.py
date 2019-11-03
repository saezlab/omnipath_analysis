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

import os
import itertools
import time

from pypath import session_mod

import omnipath2
from omnipath2 import settings as op2_settings


class PathBase(session_mod.Logger):
    
    
    timestamp = time.strftime('%Y%m%d')
    
    
    def __init__(
            self,
            fname = None,
            fname_param = (),
            fname_timestamp = True,
            timestamp_override = None,
            filetype = 'pdf',
            target_dir = None,
            dir_timestamp = True,
            **kwargs
        ):
        
        if not hasattr(self, '_logger'):
            
            session_mod.Logger.__init__(self, name = log_label or 'op2.path')
        
        for attr, val in itertools.chain(locals().items(), kwargs.items()):
            
            setattr(self, attr, val)
    
    
    def main(self):
        
        self.set_timestamp()
        self.set_directory()
        self.set_format()
        self.set_path()
    
    
    def set_timestamp(self, timestamp = None, strftime = ''):
        
        self.timestamp = (
            timestamp
                or
            time.strftime(strftime)
                or
            self.timestamp_override
                or
            omnipath2.data.timestamp
                or
            self.timestamp
        )
    
    
    def set_directory(self, target_dir = None):
        
        target_dir = target_dir or self.target_dir
        
        if hasattr(omnipath2.data, target_dir):
            
            self.target_dir = getattr(omnipath2.data, target_dir)
        
        else:
            
            self.target_dir = target_dir or op2_settings.get(target_dir)
            
            if self.dir_timestamp:
                
                self.target_dir = os.path.join(
                    self.target_dir,
                    self.timestamp,
                )
        
        os.makedirs(self.target_dir, exist_ok = True)
    
    
    def set_format(self, filetype = None):
        
        self.filetype = (
            filetype
                if filetype else
            self.filetype
                if hasattr(self, 'filetype') else
            'pdf'
        )
    
    
    def set_path(self, fname = None):
        
        self.fname = (
            fname
                or
            op2_settings.get(self.fname)
                or
            self.fname
        )
        
        if not self.fname:
            
            self._log('No output file name provided.')
            return
        
        self.fname = self.fname % self.fname_param
        self.fname = '%s%s.%s' % (
            self.fname,
            '__%s' % self.timestamp if self.fname_timestamp else '',
            self.filetype,
        )
        self.path = os.path.join(self.target_dir, self.fname)
        
        self._log('Plotting to `%s`.' % self.path)
