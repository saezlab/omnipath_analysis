#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

from pypath import websrvtab
from pypath import settings

settings.setup(
    cosmic_credentials = {
        'user': 'denes.turei@embl.de',
        'passwd': 'San33Ger##',
    },
    network_expand_complexes = False,
)

# if you want to use a non-default cache directory:
# settings.setup(cachedir = '~/.pypath/cache')

builder = websrvtab.WebserviceTables()
builder.main()
