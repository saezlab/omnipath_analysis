#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os

import matplotlib.pyplot as plt
import pandas as pd

import pypath
from pypath import intercell

# Colors!
cg = (87/255, 171/255, 39/255)
cb = (0/255, 84/255, 159/255)

cachedir = '~/pypath_cache'

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

pypath.settings.setup(cachedir=cachedir)

i = intercell.IntercellAnnotation()

dir(i)

len(i.receptors)
[(k, len(v)) for k, v in i.receptors_by_resource.items()]
