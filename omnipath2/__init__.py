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

import os
import importlib as imp
import copy

from pypath.share import settings as pp_settings
from pypath.share import session as session_mod
from omnipath2 import settings as settings_mod
from omnipath2 import database as _database_mod

if '_colors_mod' not in globals():

    from omnipath2 import colors as _colors_mod

else:

    imp.reload(_colors_mod)

if '_files_mod' not in globals():

    from omnipath2 import files as _files_mod

else:

    imp.reload(_files_mod)


if '_logger' not in globals():

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

    globals()['data'] = _database_mod.Database(**param)
    globals()['colors'] = _colors_mod.Colors()
    globals()['files'] = _files_mod.Files()


init()
