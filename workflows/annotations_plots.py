#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de

import os

import matplotlib.pyplot as plt
import pandas as pd

import pypath
from pypath import annot

# Colors!
cg = (87/255, 171/255, 39/255)
cb = (0/255, 84/255, 159/255)

if os.getcwd().endswith('omnipath2'):
    os.chdir('workflows')

if not os.path.exists('new_cache'):
    os.makedirs('new_cache')

pypath.settings.setup(cachedir='new_cache')

a = annot.AnnotationTable(keep_annotators=True)
a.load()

a.to_dataframe()

# XXX: Pending fix from pypath #75
#a.names[1]
#for i in a.names:
#    try:
#        print('__'.join(i))
#    except:
#        print('#######################')
#        print(i)
#        print('#######################')

#['__'.join(name) for name in a.names]
#dir(a)
#pd.DataFrame(a.data, index=a.uniprots, columns=a.names) ==
