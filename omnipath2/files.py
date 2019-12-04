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

import os
import json

from pypath import common
from pypath import session_mod

from omnipath2 import settings


class Files(session_mod.Logger):
    
    
    def __init__(self, json_file = None):
        
        session_mod.Logger.__init__(self, name = 'op2.files')
        
        self.json_file = json_file or settings.get('files_json')
        
        self.init()
    
    
    def init(self):
        
        self.files = {
            'recent': {},
            'history': {},
        }
        self.read_files_db()
    
    
    def read_files_db(self):
        
        if os.path.exists(self.json_file):
            
            with open(self.json_file, 'r') as fp:
                
                self.files = json.load(fp)
    
    
    def write_files_db(self):
        
        json_dump = json.dumps(
            self.files,
            sort_keys = True,
            indent = 4,
            separators = (',', ': '),
        )
        
        with open(self.json_file, 'w') as fp:
            
            fp.write(json_dump)
    
    
    def update_record(self, path):
        
        fname = os.path.basename(os.path.splitext(path)[0])
        fname = fname.split('__')[0]
        
        self._update_record(fname, path)
    
    
    def _update_record(self, fname, path):
        
        self.read_files_db()
        
        if fname in self.files['recent']:
            
            if fname not in self.files['history']:
                
                self.files['history'][fname] = set()
            
            self.files['history'][fname] = common.to_set(
                self.files['history'][fname]
            )
            self.files['history'][fname].add(
                self.files['recent'][fname][0]
                    if isinstance(self.files['recent'][fname], list) else
                self.files['recent'][fname]
            )
            self.files['history'][fname] = list(self.files['history'][fname])
        
        self.files['recent'][fname] = path
        
        self.write_files_db()
