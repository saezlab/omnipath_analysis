#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

import imp
import pprint
import collections
import iteritems

from pypath import mapping
from pypath import dataio
from pypath import intercell_annot
from pypath import annot
from pypath import intercell
from pypath import main
from pypath import network


def reload():
    
    imp.reload(intercell_annot)
    imp.reload(annot)
    imp.reload(intercell)


# this attempts to come to a consensus about the roles of the proteins
# then the next class will check these against the network resources
# and provide reports for plotting
i = intercell.IntercellAnnotation()

pa = main.PyPath()
pa.init_network(pfile = '../../omnipath_webservice/omnipath.pickle')
n = network.Network.from_igraph(pa)

pprint.pprint([(k, len(v)) for k, v in i.classes.items()])

i.export(
    fname = 'intercell_classes_20190411.tsv',
    sep = '\t',
    index = False,
)





cov_hdr = ['typ', 'cls', 'total', 'omnipath']
coverages = []
for typ, ccls in intercell_annot.class_types.items():
    
    for cls in ccls:
        
        total = len(i.classes[cls])
        in_network = len(i.classes[cls] & set(pa.graph.vs['name']))
        
        coverages.append([
            typ,
            cls,
            total,
            in_network,
        ])


with open('main_coverage.tsv', 'w') as fp:
    _ = fp.write('\t'.join(cov_hdr))
    _ = fp.write('\n')
    for l in coverages:
        _ = fp.write('%s\t%s\t%u\t%u\n' % tuple(l))


class_names = set(itertools.chain(*intercell_annot.class_types.values()))
class_types = dict(
    (cls, typ)
    for typ, ccls in intercell_annot.class_types.items()
    for cls in ccls
)
cov_hdr = ['typ', 'mainclass', 'cls', 'total', 'omnipath']
coverages_all = []
annot_resource_overlaps = collections.defaultdict(dict)

for cls, cls_proteins in i.classes.items():
    cls_split = cls.split('_')
    mainclass = None
    for j in range(len(cls_split)):
        this_part = '_'.join(cls_split[:j])
        if this_part in class_names:
            mainclass = this_part
    coverages_all.append([
        (
            class_types[mainclass]
                if mainclass in class_types else
            class_types[cls]
                if cls in class_types else
            ''
        ),
        mainclass or '',
        cls,
        len(cls_proteins),
        len(cls_proteins & set(pa.graph.vs['name'])),
    ])


# intercell network

icn = pd.merge(
    n.records,
    i.df,
    suffixes = ['', '_a'],
    how = 'inner',
    left_on = 'id_a',
    right_on = 'uniprot',
)

icn = pd.merge(
    icn,
    i.df,
    suffixes = ['_a', '_b'],
    how = 'inner',
    left_on = 'id_b',
    right_on = 'uniprot',
)


conn_hdr = ['cat0', 'cat1', 'size0', 'size1', 'conn']
connections = []

for c0, c1 in itertools.combinations_with_replacement(i.class_names, 2):
    if i.class_types[c0] == 'above_main' or i.class_types[c1] == 'above_main':
        continue
    ncon = sum(
        np.logical_or(
            np.logical_and(
                icn.category_a == c0,
                icn.category_b == c1,
            ),
            np.logical_and(
                icn.category_a == c1,
                icn.category_b == c0,
            )
        )
    )
    connections.append([
        c0,
        c1,
        len(i.classes[c0]),
        len(i.classes[c1]),
        ncon,
    ])

with open('connections.tsv', 'w') as fp:
    _ = fp.write('\t'.join(conn_hdr))
    _ = fp.write('\n')
    for conn in connections:
        _= fp.write('%s\t%s\t%u\t%u\t%u\n' % tuple(conn))

# overlaps between categories
cat_overlap_hdr = ['cat0', 'cat1', 'size0', 'size1', 'total', 'overlap']
category_overlaps = []

for c0, c1 in itertools.product(i.class_names, i.class_names):
    if i.class_types[c0] == 'above_main' or i.class_types[c1] == 'above_main':
        continue
    category_overlaps.append([
        c0,
        c1,
        len(i.classes[c0]),
        len(i.classes[c1]),
        len(i.classes[c0] | i.classes[c1]),
        len(i.classes[c0] & i.classes[c1]),
    ])

with open('category_overlaps.tsv', 'w') as fp:
    _ = fp.write('\t'.join(cat_overlap_hdr))
    _ = fp.write('\n')
    for c_ol in category_overlaps:
        _ = fp.write('%s\t%s\t%u\t%u\t%u\t%u\n' % tuple(c_ol))
