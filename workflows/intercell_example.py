#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

import imp

from pypath import mapping
from pypath import dataio
from pypath import annot
from pypath import intercell


def reload():
    
    imp.reload(annot)
    imp.reload(intercell)


# this attempts to come to a consensus about the roles of the proteins
# then the next class will check these against the network resources
# and provide reports for plotting
i = intercell.IntercellAnnotation()

len(i.receptors)
[(k, len(v)) for k, v in i.receptors_by_resource.items()]
