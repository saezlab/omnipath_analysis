#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

from pypath import websrvtab
from pypath import settings

# create a new cache directory, this you can
# later copy over the default cache `~/.pypath/cache`
settings.setup(cachedir = 'new_cache')

# call any method which downloads something
# or just run the webservice tables build process
# as this downloads almost everything:
builder = websrvtab.WebserviceTables()
builder.main()
