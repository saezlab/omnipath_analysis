#!/usr/bin/env python
#-*- coding: utf-8 -*-

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

import numpy as np
import pandas as pd

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
        self.complex_args = complex_args or {}
        self.annotation_args = annotation_args or {}
        self.intercell_args = intercell_args or {}
        
        self.omnipath_pickle = omnipath_pickle
        self.complex_pickle = complex_pickle
        self.annotation_pickle = annotation_pickle
        self.intercell_pickle = intercell_pickle
        
        self.datadir = datadir
        self.date = time.strftime('%Y%m%d')
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
        
        for attr in (
            'complex',
            'network',
            'igraph_network',
            'annot',
            'intercell',
        ):
            
            if hasattr(self, attr):
                
                getattr(self, attr).reload()
    
    
    def main(self):
        
        self.setup()
        self.load()
        self.export_tables()
    
    
    def setup(self):
        
        self.datadir = self.datadir or (
            'data' if os.path.exists('data') else os.path.join('..', 'data')
        )
    
    
    def load(self):
        
        self.load_complex()
        self.load_annot()
        self.load_intercell()
        self.load_network()
        self.build_intercell_network()
    
    
    def export_tables(self):
        
        self.export_intercell_classes()
        self.collect_classes()
        self.export_intercell_coverages()
        self.export_intercell_coverages_by_resource()
        self.build_intercell_network()
        self.export_connections()
        self.export_overlaps()
    
    
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
        
        if self.omnipath_pickle:
            
            self.igraph_network.init_network(pfile = self.omnipath_pickle)
            
        else:
            
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
                
                total = len(self.intercell.classes[cls])
                in_network = len(
                    self.intercell.classes[cls] &
                    set(self.igraph_network.graph.vs['name'])
                )
                
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
    
    
    def export_stats_by_resource(self):
        
        def add_stats_record(
                stats,
                cls0,
                cls1,
                elements0,
                elements1,
                reference_set,
                typ,
            ):
        
            parent0 = self.intercell.parents[cls0]
            parent1 = self.intercell.parents[cls1]
            
            cls0_elements = elements0 & reference_set
            cls1_elements = elements1 & reference_set
            
            overlap = cls0_elements & cls1_elements
            
            parent0_elements = classes[parent0] if parent0 else set()
            parent1_elements = classes[parent1] if parent1 else set()
            
            parent0_elements = parent0_elements & reference_set
            parent1_elements = parent1_elements & reference_set
            
            omnipath0 = cls0_elements & omnipath
            omnipath1 = cls1_elements & omnipath
            
            degree0 = np.mean([
                v.degree()
                for v in self.igraph_network.vs
                if v['name'] in cls0_elements]
            ) if typ == 'protein' else np.nan
            degree1 = np.mean([
                v.degree()
                for v in self.igraph_network.vs
                if v['name'] in cls1_elements]
            ) if typ == 'protein' else np.nan
            
            stats.append([
                cls0,
                cls1,
                labels[cls0] if cls0 in labels else '',
                labels[cls1] if cls1 in labels else '',
                class_types[cls0] if cls0 in class_types else 'sub',
                class_types[cls1] if cls1 in class_types else 'sub',
                typ,
                len(reference_set),
                len(omnipath) if typ == 'protein' else 0,
                parent0,
                parent1,
                len(cls0_elements),
                len(cls1_elements),
                len(overlap),
                len(omnipath0),
                len(omnipath1),
                len(parent0_elements),
                len(parent1_elements),
                degree0,
                degree1,
            ] + [np.nan] * 9)
        
        hdr = [
            'name_cls0',
            'name_cls1',
            'label0',
            'label1',
            'typ_cls0',
            'typ_cls1',
            'entity',
            'total',
            'omnipath',
            'parent0',
            'parent1',
            'size_cls0',
            'size_cls1',
            'overlap',
            'omnipath0',
            'omnipath1',
            'size_parent0',
            'size_parent1',
            'degree0',
            'degree1',
            'con_all',
            'con_0to1',
            'con_0to1_stim',
            'con_0to1_inh',
            'con_1to0',
            'con_1to0_stim',
            'con_1to0_inh',
            'con_in0',
            'con_in1',
        ]
        
        stats = []
        
        path = os.path.join(
            self.datadir,
            'stats_by_resource_%s.tsv' % self.date,
        )
        
        proteins = set(self.annot.proteins)
        complexes = set(self.annot.complexes)
        classes = self.intercell.classes
        labels = self.intercell.labels
        omnipath = set(self.igraph_network.vs['name'])
        class_types = self.intercell.class_types
        
        for (
            (cls0, cls0_elements),
            (cls1, cls1_elements),
        ) in (
            itertools.product(
                classes.items(),
                classes.items(),
            )
        ):
            
            add_stats_record(
                stats,
                cls0,
                cls1,
                cls0_elements,
                cls1_elements,
                proteins,
                'protein',
            )
            
            add_stats_record(
                stats,
                cls0,
                cls1,
                cls0_elements,
                cls1_elements,
                complexes,
                'complex',
            )
        
        self.stats = pd.DataFrame(
            stats,
            columns = hdr,
        )
        
        self.stats.to_csv(path, index = False, sep = '\t')
    
    
    def build_intercell_network(self):
        
        self.intercell_network = pd.merge(
            self.network.records,
            self.intercell.df,
            suffixes = ['', '_a'],
            how = 'inner',
            left_on = 'id_a',
            right_on = 'uniprot',
        )
        self.intercell_network.id_a = (
            self.intercell_network.id_a.astype('category')
        )
        
        self.intercell_network = pd.merge(
            self.intercell_network,
            self.intercell.df,
            suffixes = ['_a', '_b'],
            how = 'inner',
            left_on = 'id_b',
            right_on = 'uniprot',
        )
        self.intercell_network.id_b = (
            self.intercell_network.id_b.astype('category')
        )
    
    
    def connections_between_categories(
            self,
            cat_a,
            cat_b,
            directed = True,
            effect = None,
        ):
        
        
        idx = np.logical_and(
            np.logical_and(
                self.intercell_network.category_a == cat_a,
                self.intercell_network.category_b == cat_b,
            ),
            self.intercell_network.directed == directed,
        )
        
        if not directed:
            
            idx = np.logical_or(
                idx,
                np.logical_and(
                    np.logical_and(
                        self.intercell_network.category_a == cat_b,
                        self.intercell_network.category_b == cat_a,
                    ),
                    self.intercell_network.directed == 0,
                )
            )
            
        elif effect is not None:
            
            idx = np.logical_and(
                idx,
                self.intercell_network.effect == effect,
            )
        
        return self.intercell_network[idx]
    
    
    def count_connections_between_categories(
            self,
            cat_a,
            cat_b,
            directed = True,
            effect = None,
        ):
        
        return self.connections_between_categories(
            cat_a = cat_a,
            cat_b = cat_b,
            directed = directed,
            effect = effect,
        ).shape[0]
    
    
    def all_connections_between_categories(self, cat_a, cat_b):
        
        return (
            self.count_connections_between_categories(cat_a, cat_b) +
            self.count_connections_between_categories(cat_b, cat_a) +
            self.count_connections_between_categories(
                cat_a = cat_a,
                cat_b = cat_b,
                directed = False,
            )
        )
    
    
    def connections_a_to_b(self, cat_a, cat_b):
        
        return self.count_connections_between_categories(cat_a, cat_b)
    
    
    def a_stimulates_b(self, cat_a, cat_b):
        
        return self.count_connections_between_categories(
            cat_a = cat_a,
            cat_b = cat_b,
            effect = 1,
        )
    
    
    def a_inhibits_b(self, cat_a, cat_b):
        
        return self.count_connections_between_categories(
            cat_a = cat_a,
            cat_b = cat_b,
            effect = -1,
        )
    
    
    def add_intercell_network_stats(self):
        
        self.intercell_network.groupby(by = 'pair')
    
    
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
                self.intercell.class_types[c0] == 'above_main' or
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
                len(self.intercell.classes[c0]),
                len(self.intercell.classes[c1]),
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
        cat_overlap_hdr = [
            'cat0', 'cat1', 'size0', 'size1', 'total', 'overlap'
        ]
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
