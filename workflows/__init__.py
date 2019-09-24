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

import collections

from pypath import settings as pp_settings
from workflows import settings as settings_mod


def setup(environment):
    """
    environment: either nico or denes
    """
    
    settings_mod.setup(**getattr(settings_mod, environment))
    
    pp_cachedir = pp_settings.get('cachedir')
    
    if pp_cachedir:
        
        pp_settings.setup(cachedir = pp_cachedir)
