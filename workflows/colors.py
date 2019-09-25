#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019 Saez Lab
#
# OmniPath2 analysis and figures suit
#
# Authors:
#
# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de
#
# Dénes Türei
# turei.denes@gmail.com
#

import os
import re
import collections

from pypath import session_mod

from workflows import settings as op2_settings


recolor = re.compile(
    r'([0-9]+)\t'
    r'([0-9]+)\t'
    r'([0-9]+)\t'
    r'(\w+)'
)

class Colors(session_mod.Logger):
    
    
    def __init__(self):
        
        session_mod.Logger.__init__(self, name = 'colors')
        
        self.palette_dir = op2_settings.get('palette_dir')
        self.palettes = {}
    
    
    def read_palettes(self):
        
        for fname in os.listdir(self.palette_dir):
            
            if os.path.splitext(fname)[-1] == '.gpl':
                
                self.read_palette(os.path.join(self.palette_dir, fname))
    
    
    def read_palette(self, path, name = None):
        
        self._log('Loading color palette from `%s`.' % path)
        
        colors = collections.OrderedDict()
        
        with open(path, 'r') as fp:
            
            for l in fp:
                
                if not name and l.startswith('Name:'):
                    
                    name = l.split(':')[-1].strip()
                
                m = recolor.search(l)
                
                if not m:
                    
                    continue
                
                r, g, b, label = m.groups(0)
                
                colors[label] = tuple(int(val) / 255 for val in (r, g, b))
        
        if not name:
            
            name = os.path.basename(os.path.splitext(path)[0])
        
        self.palettes[name] = colors
        
        self._log(
            'Color palette `%s` with %u colors from `%s`.' % (
                name,
                len(colors),
                path,
            )
        )
