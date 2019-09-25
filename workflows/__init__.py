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

from pypath import settings as pp_settings
from workflows import settings as settings_mod
from workflows import compile_data
from workflows import colors


def setup(environment):
    """
    environment: either nico or denes
    """
    
    settings_mod.setup(**getattr(settings_mod, environment))
    
    pp_cachedir = pp_settings.get('cachedir')
    
    if pp_cachedir:
        
        pp_settings.setup(cachedir = pp_cachedir)


def init(environment = None, **kwargs):
    
    param = (
        copy.deepcopy(globals()['OP2_DB_ARGS'])
            if 'OP2_DB_ARGS' in globals() else
        {}
    )
    
    param.update(kwargs)
    
    environment = environment or os.path.split(os.path.expanduser('~'))[-1]
    
    setup(environment)
    
    globals()['data'] = compile_data.Database(**param)
    globals()['colors'] = colors.Colors()


init()
