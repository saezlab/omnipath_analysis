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
import json

from pypath import session_mod

from omnipath2 import settings


class Files(session_mod.Logger):
    
    
    def __init__(self, json_file = None):
        
        session_mod.Logger.__init__(self, name = 'op2.files')
        
        self.json_file = json_file or settings.get('files_json')
        
        
