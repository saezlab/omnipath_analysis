#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

from pypath import websrvtab
from pypath import settings

# if you want to use a non-default cache directory:
# settings.setup(cachedir = '~/.pypath/cache')

builder = websrvtab.WebserviceTables()
builder.main()
