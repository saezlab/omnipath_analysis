#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

from pypath import websrvtab
from pypath import settings
from pypath import annot
from pypath import complex

settings.setup(
    cosmic_credentials = {
        'user': 'denes.turei@embl.de',
        'passwd': 'qwDFbnm76#',
    },
    network_expand_complexes = False,
)

# if you want to use a non-default cache directory:
# settings.setup(cachedir = '~/.pypath/cache')


complex.init_db(pickle_file = 'complexes.pickle')
annot.init_db(pickle_file = 'annotations2.pickle')

builder = websrvtab.WebserviceTables(
    annot_args = {'pickle_file': 'annotations2.pickle'},
    complex_args = {'pickle_file': 'complexes.pickle'},
)
builder.main()
