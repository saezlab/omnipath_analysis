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
from pypath import session_mod
from omnipath2 import settings as settings_mod
from omnipath2 import database
from omnipath2 import colors


_logger = session_mod.Logger(name = 'op2.init')
_log = _logger._log


_log('Welcome to the OmniPath2 analysis and figures suite.')


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
    
    _log(
        'You can customize the database building process by '
        'setting parameters in the `OP2_DB_ARGS` global variable '
        'or by calling `init` again with keyword arguments or after '
        'setting values in the `settings` module.'
    )
    
    param.update(kwargs)
    
    environment = environment or os.path.split(os.path.expanduser('~'))[-1]
    
    setup(environment)
    
    globals()['data'] = database.Database(**param)
    globals()['colors'] = colors.Colors()


init()
