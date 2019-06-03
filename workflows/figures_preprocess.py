#!/usr/bin/env python

# Denes Turei 2019
# turei.denes@gmail.com

import imp
import os
import sys
import pprint
import copy
import collections
import itertools
import time

from pypath import mapping
from pypath import dataio
from pypath import intercell_annot
from pypath import annot
from pypath import intercell
from pypath import complex
from pypath import main
from pypath import network
from pypath import session_mod
from pypath import settings


def reload():
    
    imp.reload(intercell_annot)
    imp.reload(annot)
    imp.reload(intercell)


class FiguresPreprocess(session_mod.Logger):
    
    omnipath_args_default = {
        'kinase_substrate_extra': True,
        'ligand_receptor_extra': True,
    }
    
    def __init__(
            self,
            omnipath_pickle = None,
            complex_pickle = None,
            annotation_pickle = None,
            intercell_pickle = None,
            omnipath_args = None,
            complex_args = None,
            annotation_args = None,
            intercell_args = None,
            datadir = None,
        ):
        
        session_mod.Logger.__init__(self, name = 'figures_preproc')
        
        self.omnipath_args = copy.deepcopy(self.omnipath_args_default)
        self.omnipath_args.update(omnipath_args or {})
        
        if omnipath_pickle and os.path.exists(omnipath_pickle):
            
            self.omnipath_args['pfile'] = omnipath_pickle
        
        self.complex_args = complex_args or {}
        self.annotation_args = annotation_args or {}
        self.intercell_args = intercell_args or {}
        self.complex_pickle = complex_pickle
        self.annotation_pickle = annotation_pickle
        self.intercell_pickle = intercell_pickle
        
        self.datadir = datadir
        self.date = time.strftime('%Y%m%d')
    
    
    def main(self):
        
        self.setup()
        self.load()
    
    
    def setup(self):
        
        self.datadir = self.datadir or (
            'data' if os.path.exists('data') else os.path.join('..', 'data')
        )
    
    
    def load(self):
        
        self.load_complex()
        self.load_annot()
        self.load_intercell()
        self.load_network()
    
    
    def load_complex(self):
        
        self.complex = complex.get_db(
            pickle_file = self.complex_pickle,
            **self.complex_args
        )
    
    
    def load_annot(self):
        
        self.annot = annot.get_db(
            pickle_file = self.annotation_pickle,
            **self.annotation_args
        )
    
    
    def load_intercell(self):
        
        self.intercell = intercell.IntercellAnnotation(
            pickle_file = self.intercell_pickle,
            **self.intercell_args
        )
    
    
    def load_network(self):
        
        self.igraph_network = main.PyPath()
        self.igraph_network.load_omnipath(**self.omnipath_args)
        self.network = network.Network.from_igraph(self.igraph_network)
    
    
    def print_intercell_classes(self):
        
        pprint.pprint([
            (k, len(v))
            for k, v in self.intercell.classes.items()
        ])
    
    
    def export_intercell_classes(self):
        
        self.intercell.export(
            fname = os.path.join(
                self.datadir,
                'intercell_classes_%s.tsv' % self.date,
            ),
            sep = '\t',
            index = False,
        )
    
    
    def export_intercell_coverages(self):
        
        cov_hdr = ['typ', 'cls', 'total', 'omnipath']
        self.intercell_coverages = []
        
        for typ, ccls in intercell_annot.class_types.items():
            
            for cls in ccls:
                
                total = len(i.classes[cls])
                in_network = len(i.classes[cls] & set(pa.graph.vs['name']))
                
                self.intercell_coverages.append([
                    typ,
                    cls,
                    total,
                    in_network,
                ])
        
        path = os.path.join(
            self.datadir,
            'main_coverage_%s.tsv' % self.date,
        )
        
        with open(path, 'w') as fp:
            
            _ = fp.write('\t'.join(cov_hdr))
            _ = fp.write('\n')
            
            for l in self.intercell_coverages:
                
                _ = fp.write('%s\t%s\t%u\t%u\n' % tuple(l))
    
    
    def collect_classes(self):
        
        self.class_names = set(
            itertools.chain(
                *intercell_annot.class_types.values()
            )
        )
        
        self.class_types = dict(
            (cls, typ)
            for typ, ccls in intercell_annot.class_types.items()
            for cls in ccls
        )
    
    
    def export_intercell_coverages2(self):
        """
        I have no idea why 2 methods
        """
        
        cov_hdr = ['typ', 'mainclass', 'cls', 'total', 'omnipath']
        self.coverages_all = []
        annot_resource_overlaps = collections.defaultdict(dict)
        
        path = os.path.join(
            self.datadir,
            'main_coverage_2_%s.tsv' % self.date,
        )
        
        for cls, cls_proteins in i.classes.items():
            
            cls_split = cls.split('_')
            mainclass = None
            
            for j in range(len(cls_split)):
                
                this_part = '_'.join(cls_split[:j])
                
                if this_part in class_names:
                    
                    mainclass = this_part
            
            self.coverages_all.append([
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
    
    
    def build_intercell_network(self):
        
        self.intercell_network = pd.merge(
            self.network.records,
            self.intercell.df,
            suffixes = ['', '_a'],
            how = 'inner',
            left_on = 'id_a',
            right_on = 'uniprot',
        )
        
        self.intercell_network = pd.merge(
            self.intercell_network,
            self.intercell.df,
            suffixes = ['_a', '_b'],
            how = 'inner',
            left_on = 'id_b',
            right_on = 'uniprot',
        )
    
    
    def export_connections(self):
        
        conn_hdr = ['cat0', 'cat1', 'size0', 'size1', 'conn']
        self.connections = []
        
        for c0, c1 in (
            itertools.combinations_with_replacement(
                self.intercell.class_names,
                2,
            )
        ):
            
            if (
                self.intercell.i.class_types[c0] == 'above_main' or
                self.intercell.class_types[c1] == 'above_main'
            ):
                continue
            
            numof_connections = sum(
                np.logical_or(
                    np.logical_and(
                        self.intercell_network.category_a == c0,
                        self.intercell_network.category_b == c1,
                    ),
                    np.logical_and(
                        self.intercell_network.category_a == c1,
                        self.intercell_network.category_b == c0,
                    )
                )
            )
            
            self.connections.append([
                c0,
                c1,
                len(i.classes[c0]),
                len(i.classes[c1]),
                numof_connections,
            ])
        
        path = os.path.join(
            self.datadir,
            'connections_%s.tsv' % self.date,
        )
        
        with open(path, 'w') as fp:
            
            _ = fp.write('\t'.join(conn_hdr))
            _ = fp.write('\n')
            
            for conn in self.connections:
                
                _= fp.write('%s\t%s\t%u\t%u\t%u\n' % tuple(conn))
    
    
    def export_overlaps(self):
        
        # overlaps between categories
        cat_overlap_hdr = ['cat0', 'cat1', 'size0', 'size1', 'total', 'overlap']
        self.category_overlaps = []
        
        for c0, c1 in (
            itertools.product(
                self.intercell.class_names,
                self.intercell.class_names,
            )
        ):
            
            if (
                self.intercell.class_types[c0] == 'above_main' or
                self.intercell.class_types[c1] == 'above_main'
            ):
                continue
            
            self.category_overlaps.append([
                c0,
                c1,
                len(self.intercell.classes[c0]),
                len(self.intercell.classes[c1]),
                len(self.intercell.classes[c0] | self.intercell.classes[c1]),
                len(self.intercell.classes[c0] & self.intercell.classes[c1]),
            ])
        
        path = os.path.join(
            self.datadir,
            'category_overlaps_%s.tsv' % self.date,
        )
        
        with open(path, 'w') as fp:
            
            _ = fp.write('\t'.join(cat_overlap_hdr))
            _ = fp.write('\n')
            
            for c_ol in self.category_overlaps:
                
                _ = fp.write('%s\t%s\t%u\t%u\t%u\t%u\n' % tuple(c_ol))
