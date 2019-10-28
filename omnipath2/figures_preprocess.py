#!/usr/bin/env python
#-*- coding: utf-8 -*-

#
#  This file is part of the `omnipath2` Python module
#
#  Copyright
#  2019
#  Heidelberg University, Uniklinik RWTH Aachen
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from future.utils import iteritems

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
from pypath import progress
from pypath import data_formats
from pypath import intera

import omnipath2
from omnipath2 import settings as op2_settings

def reload():
    
    imp.reload(intercell_annot)
    imp.reload(annot)
    imp.reload(intercell)


class FiguresPreprocess(session_mod.Logger):
    
    
    def __init__(self, network_dataset = 'omnipath'):
        
        session_mod.Logger.__init__(self, name = 'op2.fig_preproc')
        
        self.data = omnipath2.data
        self.date = time.strftime('%Y%m%d')
        self.network_dataset = network_dataset
    
    
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
        self.build()
        self.export_tables()
    
    
    def setup(self):
        
        self.figures_dir = self.data.figures_dir
        self.tables_dir = self.data.tables_dir
    
    
    def load(self):
        
        self.load_complex()
        self.load_annot()
        self.load_intercell()
        self.load_network()
    
    
    def build(self):
        
        self.count_connections()
        self.build_intercell_network()
        self.count_connections_groupwise()
        self.count_connections_pairwise()
    
    
    def export_tables(self):
        
        self.export_intercell_classes()
        self.collect_classes()
        self.export_intercell_coverages()
        self.export_intercell_coverages_by_resource()
        self.export_connections()
        self.export_overlaps()
    
    
    def load_complex(self):
        
        self._log('Loading the complex database.')
        self.complex = self.data.get_db('complex')
    
    
    def load_annot(self):
        
        self._log('Loading the annotation database.')
        self.annot = self.data.get_db('annotations')
    
    
    def load_intercell(self):
        
        self._log('Loading the intercell annotation database.')
        self.intercell = self.data.get_db('intercell')
    
    
    def load_network(self):
        
        self._log('Loading the signaling network.')
        self.igraph_network = self.data.get_db(self.network_dataset)
        self.network = self.data.network_df(self.network_dataset)
        self.network_by_source = (
            self.data.network_df_by_source(self.network_dataset)
        )
    
    
    def set_network(self, dataset, by_source = False):
        
        self.data.set_network(dataset = dataset, by_source = by_source)
        self.network = self.data.intercell.network
    
    
    def count_connections(self):
        
        self._log('Counting connections in the network.')
        
        self.con_omnipath = len(
            set(
                zip(
                    self.network.records.id_a,
                    self.network.records.id_b,
                )
            )
        )
        self.con_omnipath_undir = (
            np.logical_not(self.network.records.directed).sum()
        )
        directed = self.network.records[self.network.records.directed]
        self.con_omnipath_dir = len(set(zip(directed.id_a, directed.id_b)))
        self.con_omnipath_stim = (self.network.records.effect == 1).sum()
        self.con_omnipath_inh = (self.network.records.effect == -1).sum()
    
    
    def print_intercell_classes(self):
        
        pprint.pprint([
            (k, len(v))
            for k, v in self.intercell.classes.items()
        ])
    
    
    def export_intercell_classes(self):
        
        self.intercell.export(
            fname = os.path.join(
                self.tables_dir,
                op2_settings.get('intercell_classes_tsv') % self.date,
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
            self.tables_dir,
            op2_settings.get('main_coverage_tsv') % self.date,
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
    
    
    def count_connections_groupwise(self):
        """
        Counts the overall number of connections for each category.
        """
        
        def count(df, undirected = False, in_out = 'IN'):
            
            sorter = sorted if undirected else lambda x: x
            
            query = (
                'category_a == "{}" | category_b == "{}"'
                    if undirected else
                'category_a == "{}"'
                    if in_out == 'OUT' else
                'category_b == "{}"'
            )
            
            return collections.defaultdict(
                int,
                dict(
                    (
                        cls,
                        len(
                            set(
                                map(
                                    tuple,
                                    map(
                                        sorter,
                                        zip(
                                            df.query(
                                                query.format(
                                                    *(cls,) *
                                                    query.count('{}')
                                                )
                                            ).id_a,
                                            df.query(
                                                query.format(
                                                    *(cls,) *
                                                    query.count('{}')
                                                )
                                            ).id_b,
                                        )
                                    )
                                )
                            )
                        )
                    )
                    for cls in self.intercell.class_names
                )
            )
        
        
        self.degree_total = count(
            self.intercell_network,
            undirected = True,
        )
        
        self.degree_undirected = count(
            self.intercell_network[
                np.logical_not(self.intercell_network.directed)
            ],
            undirected = True,
        )
        
        for mode in ('IN', 'OUT'):
            
            setattr(
                self,
                'degree_%s' % mode.lower(),
                count(
                    self.intercell_network[
                        self.intercell_network.directed
                    ],
                    in_out = mode,
                ),
            )
            
            setattr(
                self,
                'degree_%s_stim' % mode.lower(),
                count(
                    self.intercell_network[
                        self.intercell_network.effect == 1
                    ],
                    in_out = mode,
                ),
            )
            
            setattr(
                self,
                'degree_%s_inh' % mode.lower(),
                count(
                    self.intercell_network[
                        self.intercell_network.effect == -1
                    ],
                    in_out = mode,
                ),
            )
    
    
    def count_connections_pairwise_old(self):
        """
        Counts the connections between categories pairwise.
        """
        
        def count(df, undirected = False):
            
            #sorter = sorted if undirected else lambda x: x
            
            con = zip(
                df.category_a,
                df.category_b
            )
            
            if undirected:
                
                con = itertools.chain(
                    con,
                    zip(df.category_b, df.category_a),
                )
            
            return collections.defaultdict(
                int,
                collections.Counter(con).items(),
            )
        
        
        self.con_directed_o    = count(
            self.intercell_network[self.intercell_network.directed]
        )
        self.con_stimulation_o = count(
            self.intercell_network[self.intercell_network.effect ==  1]
        )
        self.con_inhibition_o  = count(
            self.intercell_network[self.intercell_network.effect == -1]
        )
        self.con_undirected_o  = count(
            self.intercell_network[
                np.logical_not(self.intercell_network.directed)
            ],
            undirected = True,
        )
    
    
    def count_connections_pairwise(self):
        """
        Counts the connections between categories pairwise.
        """
        
        def counts_dict(sets_dict):
            
            return collections.defaultdict(
                int,
                dict(
                    (cats, len(elements))
                    for cats, elements in sets_dict.items()
                )
            )
        
        
        self._log('Counting connections between classes.')
        
        con_all = collections.defaultdict(set)
        undirected = collections.defaultdict(set)
        directed = collections.defaultdict(set)
        stimulation = collections.defaultdict(set)
        inhibition = collections.defaultdict(set)
        
        prg = progress.Progress(
            self.intercell_network.shape[0],
            'Counting connections',
            1000,
        )
        
        for i in self.intercell_network.itertuples():
            
            prg.step()
            
            sorted_tuple = tuple(sorted((i.id_a, i.id_b)))
            
            if not i.directed:
                
                undirected[(i.category_a, i.category_b)].add(sorted_tuple)
                undirected[(i.category_b, i.category_a)].add(sorted_tuple)
                
            else:
                
                directed[(i.category_a, i.category_b)].add((i.id_a, i.id_b))
                
                if i.effect == 1:
                    
                    stimulation[(i.category_a, i.category_b)].add(
                        (i.id_a, i.id_b)
                    )
                    
                elif i.effect == -1:
                    
                    inhibition[(i.category_a, i.category_b)].add(
                        (i.id_a, i.id_b)
                    )
            
            con_all[(i.category_a, i.category_b)].add(sorted_tuple)
            con_all[(i.category_b, i.category_a)].add(sorted_tuple)
            
        
        prg.terminate()
        
        self.con_undirected = counts_dict(undirected)
        self.con_directed = counts_dict(directed)
        self.con_stimulation = counts_dict(stimulation)
        self.con_inhibition = counts_dict(inhibition)
        self.con_all = counts_dict(con_all)
        
        self._log('Finished counting connections between classes.')
    
    
    def export_stats_by_category_pairs(self):
        """
        Exports a number of statistics for all pairs of categories.
        """
        
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
            
            cls01 = (cls0, cls1)
            cls10 = (cls1, cls0)
            
            stats.append(
                ClassesPairwiseRecord(
                    name_cls0 = cls0,
                    name_cls1 = cls1,
                    cls_label0 = (
                        class_labels[cls0] if cls0 in class_labels else ''
                    ),
                    cls_label1 = (
                        class_labels[cls1] if cls1 in class_labels else ''
                    ),
                    src_label0 = (
                        resource_labels[cls0]
                            if cls0 in resource_labels else
                        ''
                    ),
                    src_label1 = (
                        resource_labels[cls1]
                            if cls1 in resource_labels else
                        ''
                    ),
                    typ_cls0 = (
                        class_types[cls0] if cls0 in class_types else 'sub'
                    ),
                    typ_cls1 = (
                        class_types[cls1] if cls1 in class_types else 'sub'
                    ),
                    entity = typ,
                    total = len(reference_set), # total number of all proteins
                                                # or all complexes
                    omnipath = len(omnipath) if typ == 'protein' else 0,
                    parent0 = parent0,
                    parent1 = parent1,
                    # sizes
                    size_cls0 = len(cls0_elements),
                    size_cls1 = len(cls1_elements),
                    overlap_cls01 = len(overlap),
                    in_omnipath_cls0 = len(omnipath0),
                    in_omnipath_cls1 = len(omnipath1),
                    size_parent0 = len(parent0_elements),
                    size_parent1 = len(parent1_elements),
                    # connections
                    con_all = self.con_all[cls01],
                    con_0to1 = self.con_directed[cls01],
                    con_0to1_stim = self.con_stimulation[cls01],
                    con_0to1_inh = self.con_inhibition[cls01],
                    con_1to0 = self.con_directed[cls10],
                    con_1to0_stim = self.con_stimulation[cls10],
                    con_1to0_inh = self.con_inhibition[cls10],
                    # sum degrees
                    deg_total0 = self.degree_total[cls0],
                    deg_total1 = self.degree_total[cls1],
                    deg_undir0 = self.degree_undirected[cls0],
                    deg_undir1 = self.degree_undirected[cls1],
                    deg_in0 = self.degree_in[cls0],
                    deg_in1 = self.degree_in[cls1],
                    deg_out0 = self.degree_out[cls0],
                    deg_out1 = self.degree_out[cls1],
                    deg_in0_stim = self.degree_in_stim[cls0],
                    deg_in1_stim = self.degree_in_stim[cls1],
                    deg_out0_stim = self.degree_out_stim[cls0],
                    deg_out1_stim = self.degree_out_stim[cls1],
                    deg_in0_inh = self.degree_in_inh[cls0],
                    deg_in1_inh = self.degree_in_inh[cls1],
                    deg_out0_inh = self.degree_out_inh[cls0],
                    deg_out1_inh = self.degree_out_inh[cls1],
                    con_omnipath = self.con_omnipath,
                    con_omnipath_undir = self.con_omnipath_undir,
                    con_omnipath_dir = self.con_omnipath_dir,
                    con_omnipath_stim = self.con_omnipath_stim,
                    con_omnipath_inh = self.con_omnipath_inh,
                )
            )
        
        
        ClassesPairwiseRecord = collections.namedtuple(
            'ClassesPairwiseRecord',
            [
                'name_cls0',
                'name_cls1',
                'cls_label0',
                'cls_label1',
                'src_label0',
                'src_label1',
                'typ_cls0',
                'typ_cls1',
                'entity',
                'total',
                'omnipath',
                'parent0',
                'parent1',
                'size_cls0',
                'size_cls1',
                'overlap_cls01',
                'in_omnipath_cls0',
                'in_omnipath_cls1',
                'size_parent0',
                'size_parent1',
                'con_all',
                'con_0to1',
                'con_0to1_stim',
                'con_0to1_inh',
                'con_1to0',
                'con_1to0_stim',
                'con_1to0_inh',
                'deg_total0',
                'deg_total1',
                'deg_undir0',
                'deg_undir1',
                'deg_in0',
                'deg_in1',
                'deg_out0',
                'deg_out1',
                'deg_in0_stim',
                'deg_in1_stim',
                'deg_out0_stim',
                'deg_out1_stim',
                'deg_in0_inh',
                'deg_in1_inh',
                'deg_out0_inh',
                'deg_out1_inh',
                'con_omnipath',
                'con_omnipath_undir',
                'con_omnipath_dir',
                'con_omnipath_stim',
                'con_omnipath_inh',
            ],
        )
        
        stats = []
        
        path = os.path.join(
            self.tables_dir,
            op2_settings.get('intercell_network_by_resource_tsv') % self.date,
        )
        
        proteins = set(self.annot.proteins)
        complexes = set(self.annot.complexes)
        classes = self.intercell.classes
        class_labels = self.intercell.class_labels
        resource_labels = self.intercell.resource_labels
        omnipath = (
            set(self.network.records.id_a) |
            set(self.network.records.id_b)
        )
        class_types = self.intercell.class_types
        
        self._log('Building intercell network by resource table.')
        
        prg = progress.Progress(
            len(classes.items()) ** 2 / 2,
            'Collecting intercell categories stats',
            1,
            off = False,
        )
        
        i_ = 0
        
        for (
            (cls0, cls0_elements),
            (cls1, cls1_elements),
        ) in itertools.combinations_with_replacement(classes.items(), 2):
            
            i_ += 1
            
            prg.step()
            
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
            
            if i_ == 10:
                
                pass
                # break
        
        prg.terminate()
        
        self.stats = pd.DataFrame(
            stats,
            columns = stats[0]._fields,
        )
        
        self.stats.to_csv(path, index = False, sep = '\t')
        
        self._log('Finished intercell network by resource table.')
    
    
    def export_annotations_by_entity(self):
        
        all_proteins = set(dataio.all_uniprots(swissprot = True))
        all_complexes = set(self.complex.complexes.keys())
        annotations = list(self.intercell.classes.keys())
        resources = list(self.network.resources)
        
        AnnotationRecord = collections.namedtuple(
            'AnnotationRecord',
            [
                'entity_id',
                'is_complex',
                'cls',
                'parent',
                'resource_label',
                'class_label',
                'class_type',
            ],
        )
        
        tbl = []
        
        path = os.path.join(
            self.datadir,
            'annotations_by_entity_%s.tsv' % self.date,
        )
        
        for cls, elements in iteritems(self.intercell.classes):
            
            for elem in elements:
                
                tbl.append(
                    AnnotationRecord(
                        entity_id = elem.__str__(),
                        is_complex = isinstance(elem, intera.Complex),
                        cls = cls,
                        parent = self.intercell.parents[cls],
                        resource_label = self.get_resource_label(cls),
                        class_label = self.get_class_label(cls),
                        class_type = self.get_class_type(cls),
                    )
                )
        
        df = pd.DataFrame(tbl, columns = AnnotationRecord._fields)
        
        df.to_csv(path, index = False, sep = '\t')
        
        self.annotations_by_entity = df
    
    
    def export_resources_by_entity(self):
        
        ResourceRecord = collections.namedtuple(
            'ResourceRecord',
            [
                'entity_id',
                'resource',
                'is_complex',
                'resource_class',
            ],
        )
        
        
        all_entities = set(
            itertools.chain(
                *self.network.records[['id_a', 'id_b']].values
            )
        )
        
        entities_by_resource = self.network.entities_by_resource()
        
        tbl = []
        
        path = os.path.join(
            self.tables_dir,
            op2_settings.get('resources_by_entity_tsv') % self.date,
        )
        
        for resource, entities in iteritems(entities_by_resource):
            
            for entity_id in entities:
                
                tbl.append(
                    ResourceRecord(
                        entity_id = entity_id,
                        resource = resource,
                        is_complex = entity_id.startswith('COMPLEX'),
                        resource_class = (
                            data_formats.categories[resource]
                                if resource in data_formats.categories else
                            'z'
                        ),
                    )
                )
        
        df = pd.DataFrame(tbl, columns = ResourceRecord._fields)
        
        df.to_csv(path, index = False, sep = '\t')
        
        self.resources_by_entity = df
    
    
    def export_complexes_by_resource(self):
        
        ComplexRecord = collections.namedtuple(
            'ComplexRecord',
            [
                'complex_id',
                'resource',
            ],
        )
        
        tbl = []
        
        path = os.path.join(
            self.tables_dir,
            op2_settings.get('complexes_by_resource_tsv') % self.date,
        )
        
        for cplex_name, cplex in iteritems(self.complex.complexes):
            
            for resource in cplex.sources:
                
                tbl.append(
                    ComplexRecord(
                        complex_id = cplex.__str__(),
                        resource = resource,
                    )
                )
        
        df = pd.DataFrame(tbl, columns = ComplexRecord._fields)
        
        df.to_csv(path, index = False, sep = '\t')
        
        self.complexes_by_resource = df
    
    
    def get_class_type(self, cls):
        
        return (
            self.intercell.class_types[cls]
                if cls in self.intercell.class_types else
            'sub'
        )
    
    
    def get_resource_label(self, cls):
        
        return (
            self.intercell.resource_labels[cls]
                if cls in self.intercell.resource_labels else
            ''
        )
    
    
    def get_class_label(self, cls):
        
        return (
            self.intercell.class_labels[cls]
                if cls in self.intercell.class_labels else
            ''
        )
    
    
    def build_intercell_network(self, **kwargs):
        
        self._log('Building the intercell communication network.')
        
        self.intercell_network = self.intercell.network_df(**kwargs)
    
    
    def connections_of_category(
            self,
            category,
            directed = True,
            in_out = None,
            count = False,
        ):
        
        idx = np.array([False] * self.intercell_network.shape[0])
        
        if not directed or in_out != 'IN':
            
            idx = np.logical_or(
                idx,
                self.intercell_network.category_a == category,
            )
        
        if not directed or in_out == 'IN':
            
            idx = np.logical_or(
                idx,
                self.intercell_network.category_b == category,
            )
        
        return idx.sum() if count else self.intercell_network[idx]
    
    
    def count_connections_of_category(
            self,
            category,
            directed = True,
            in_out = None,
        ):
        
        edges = self.connections_of_category(
            category = category,
            directed = directed,
            in_out = in_out,
        )
        
        return sum(edges.groupby('id_a').nunique('id_b'))
    
    
    def category_out(self, category):
        
        return self.count_connections_of_category(
            category = category,
            in_out = 'OUT',
        )
    
    
    def category_in(self, category):
        
        return self.count_connections_of_category(
            category = category,
            in_out = 'IN',
        )
    
    
    def connections_between_categories(
            self,
            cat_a,
            cat_b,
            directed = True,
            effect = None,
            count = False,
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
        
        return idx.sum() if count else self.intercell_network[idx]
    
    
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
            count = True,
        )
    
    
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
            self.tables_dir,
            op2_settings.get('connections_tsv') % self.date,
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
            self.tables_dir,
            op2_settings.get('category_overlaps_tsv') % self.date,
        )
        
        with open(path, 'w') as fp:
            
            _ = fp.write('\t'.join(cat_overlap_hdr))
            _ = fp.write('\n')
            
            for c_ol in self.category_overlaps:
                
                _ = fp.write('%s\t%s\t%u\t%u\t%u\t%u\n' % tuple(c_ol))
