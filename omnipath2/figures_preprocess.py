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
import gc
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
from pypath import common
from pypath import entity

import omnipath2
from omnipath2 import settings as op2_settings
from omnipath2 import table

def reload():
    
    imp.reload(intercell_annot)
    imp.reload(annot)
    imp.reload(intercell)


class InterClassConnections(omnipath2.table.TableBase):
    
    
    def __init__(
            self,
            network_dataset = 'omnipath',
            mode = 'undirected',
            **kwargs
        ):
        
        self.network_dataset = network_dataset
        self.mode = mode
        
        param = {
            'fname': 'connections_tsv',
            'fname_param': (
                self.network_dataset,
                self.mode,
            ),
            'class_levels': {
                'main',
                'small_main',
            },
            'header': [
                'cat0',
                'cat1',
                'label0',
                'label1',
                'size0',
                'size1',
                'conn',
            ],
        }
        param.update(kwargs)
        
        omnipath2.table.TableBase.__init__(self, **param)
    
    
    def load(self):
        
        self._log(
            'Counting `%s` inter class connections in network `%s`.' % (
                self.mode,
                self.network_dataset,
            )
        )
        
        self.intercell = omnipath2.data.get_db('intercell')
        self.intercell.register_network(
            omnipath2.data.network_df(self.network_dataset)
        )
        
        mode = '' if self.mode == 'undirected' else '_%s' % self.mode
        method = 'count_inter_class_connections%s' % mode
        
        conn_hdr = ['cat0', 'cat1', 'size0', 'size1', 'conn']
        self.data = []
        
        iterator = (
            itertools.combinations_with_replacement(
                self.intercell.class_names,
                2,
            )
                if self.mode == 'undirected' else
            itertools.product(
                self.intercell.class_names,
                self.intercell.class_names,
            )
        )
        
        for c0, c1 in iterator:
            
            if (
                self.intercell.class_types[c0] not in self.class_levels or
                self.intercell.class_types[c1] not in self.class_levels
            ):
                continue
            
            numof_connections = getattr(
                self.intercell,
                method,
            )(
                source_classes = c0,
                target_classes = c1,
            )
            
            self.data.append([
                c0,
                c1,
                self.intercell.class_labels[c0],
                self.intercell.class_labels[c1],
                len(self.intercell.classes[c0]),
                len(self.intercell.classes[c1]),
                numof_connections,
            ])


class IntercellClasses(omnipath2.table.TableBase):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'fname': 'intercell_classes_tsv',
            'header': [
                'cat0',
                'cat1',
                'label0',
                'label1',
                'size0',
                'size1',
                'conn',
            ],
        }
        param.update(kwargs)
        
        omnipath2.table.TableBase.__init__(self, **param)
    
    
    def load(self):
        
        self.intercell = omnipath2.data.get_db('intercell')
        self.intercell.make_df()
        self.data = self.intercell.df
        self.header = self.data.columns


class IntercellCoverages(omnipath2.table.TableBase):
    
    
    def __init__(self, network_dataset = 'omnipath', **kwargs):
        
        self.network_dataset = network_dataset
        
        param = {
            'fname': 'main_coverage_tsv',
            'fname_param': (
                self.network_dataset,
            ),
            'header': [
                'typ',
                'cls',
                'cls_label',
                'resource',
                'entity_type',
                'total',
                'omnipath',
            ],
        }
        param.update(kwargs)
        
        omnipath2.table.TableBase.__init__(self, **param)
    
    
    def load(self):
        
        self.intercell = omnipath2.data.get_db('intercell')
        self.intercell.make_df()
        self.network = omnipath2.data.get_db(self.network_dataset)
        
        network_entities = {
            'protein': self.network.protein_entities(),
            'complex': self.network.complex_entities(),
        }
        
        self.data = []
        
        for cls in self.intercell.classes.keys():
            
            for entity_type in ('protein', 'complex'):
                
                members = self.intercell.get_class(
                    cls,
                    entity_type = entity_type,
                )
                total = len(members)
                in_network = len(
                    members &
                    network_entities[entity_type]
                )
                typ = (
                    self.intercell.class_types[cls]
                        if cls in self.intercell.class_types else
                    'sub'
                )
                res = (
                    self.intercell.resource_labels[cls]
                        if cls in self.intercell.resource_labels else
                    ''
                )
                
                self.data.append([
                    typ,
                    cls,
                    self.intercell.class_labels[cls],
                    res,
                    entity_type,
                    total,
                    in_network,
                ])


class IntercellNetworkCounts(omnipath2.table.TableBase):
    
    
    def __init__(
            self,
            network_dataset = 'omnipath',
            only_proteins = True,
            class_levels = None,
            **kwargs
        ):
        
        self.network_dataset = network_dataset
        self.only_proteins = only_proteins
        self.class_levels = common.to_set(class_levels)
        
        param = {
            'fname': 'intercell_network_by_resource_tsv',
            'fname_param': (
                self.network_dataset,
                'proteins' if self.only_proteins else 'all-entities',
                (
                    '-'.join(sorted(self.class_levels))
                        if self.class_levels else
                    'all-class-levels'
                ),
            ),
        }
        param.update(kwargs)
        
        omnipath2.table.TableBase.__init__(self, **param)
    
    
    def load(self):
        """
        Creates a table with  a number of statistics for all pairs of
        categories.
        """
        
        def add_stats_record(
                cls0,
                cls1,
                reference_set,
            ):
            
            entity_type = 'protein' if self.only_proteins else None
            
            cls0_elements = self.intercell.get_class(
                cls0,
                entity_type = entity_type,
            )
            cls1_elements = self.intercell.get_class(
                cls1,
                entity_type = entity_type,
            )
            
            parent0 = self.intercell.parents[cls0]
            parent1 = self.intercell.parents[cls1]
            
            overlap = cls0_elements & cls1_elements
            
            parent0_elements = (
                self.intercell.get_class(
                    parent0,
                    entity_type = entity_type,
                )
                if parent0 else
                set()
            )
            parent1_elements = (
                self.intercell.get_class(
                    parent1,
                    entity_type = entity_type,
                )
                if parent1 else
                set()
            )
            
            in_network = self.network.entities(entity_type = entity_type)
            in_network0 = in_network & cls0_elements
            in_network1 = in_network & cls1_elements
            
            cls01 = (cls0, cls1)
            cls10 = (cls1, cls0)
            
            class_labels = self.intercell.class_labels
            resource_labels = self.intercell.resource_labels
            class_types = self.intercell.class_types
            
            self.data.append(
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
                    total = len(reference_set), # total number of all proteins
                                                # or all complexes
                    network = len(in_network),
                    parent0 = parent0,
                    parent1 = parent1,
                    # sizes
                    size_cls0 = len(cls0_elements),
                    size_cls1 = len(cls1_elements),
                    overlap_cls01 = len(overlap),
                    in_network_cls0 = len(in_network0),
                    in_network_cls1 = len(in_network1),
                    size_parent0 = len(parent0_elements),
                    size_parent1 = len(parent1_elements),
                    # connections
                    con_all = self.con_all[cls01],
                    con_0to1 = self.con_directed[cls01],
                    con_0to1_stim = self.con_stimulatory[cls01],
                    con_0to1_inh = self.con_inhibitory[cls01],
                    con_1to0 = self.con_directed[cls10],
                    con_1to0_stim = self.con_stimulatory[cls10],
                    con_1to0_inh = self.con_inhibitory[cls10],
                    # sum degrees
                    deg_total0 = self.degree_total[cls0],
                    deg_total1 = self.degree_total[cls1],
                    deg_undir0 = self.degree_undirected_out[cls0],
                    deg_undir1 = self.degree_undirected_out[cls1],
                    deg_in0 = self.degree_directed_in[cls0],
                    deg_in1 = self.degree_directed_in[cls1],
                    deg_out0 = self.degree_directed_out[cls0],
                    deg_out1 = self.degree_directed_out[cls1],
                    deg_in0_stim = self.degree_stimulatory_in[cls0],
                    deg_in1_stim = self.degree_stimulatory_in[cls1],
                    deg_out0_stim = self.degree_stimulatory_out[cls0],
                    deg_out1_stim = self.degree_stimulatory_out[cls1],
                    deg_in0_inh = self.degree_inhibitory_in[cls0],
                    deg_in1_inh = self.degree_inhibitory_in[cls1],
                    deg_out0_inh = self.degree_inhibitory_out[cls0],
                    deg_out1_inh = self.degree_inhibitory_out[cls1],
                    con_network = self.con_network,
                    con_network_undir = self.con_network_undir,
                    con_network_dir = self.con_network_dir,
                    con_network_stim = self.con_network_stim,
                    con_network_inh = self.con_network_inh,
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
                'total',
                'network',
                'parent0',
                'parent1',
                'size_cls0',
                'size_cls1',
                'overlap_cls01',
                'in_network_cls0',
                'in_network_cls1',
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
                'con_network',
                'con_network_undir',
                'con_network_dir',
                'con_network_stim',
                'con_network_inh',
            ],
        )
        
        self.setup_data()
        
        self.count_connections_pairwise()
        self.count_connections_groupwise()
        self.count_connections()
        self.intercell.unset_interclass_network_df()
        gc.collect()
        
        self.data = []
        
        annot_entities = (
            self.annot.proteins
                if self.only_proteins else
            self.annot.reference_set
        )
        classes = self.intercell.classes
        class_labels = self.intercell.class_labels
        resource_labels = self.intercell.resource_labels
        network = (
            set(self.network_df.records.id_a) |
            set(self.network_df.records.id_b)
        )
        class_types = self.intercell.class_types
        
        self._log('Building the intercell network by resource table.')
        
        prg = progress.Progress(
            len(classes.items()) ** 2 / 2,
            'Collecting intercell categories stats',
            1,
            off = False,
        )
        
        i_ = 0
        
        for cls0, cls1 in (
            itertools.combinations_with_replacement(classes.keys(), 2)
        ):
            
            i_ += 1
            prg.step()
            
            add_stats_record(
                cls0,
                cls1,
                annot_entities,
            )
            
            if i_ == 10:
                
                pass
                # break
        
        prg.terminate()
        
        self.data = pd.DataFrame(
            self.data,
            columns = self.data[0]._fields,
        )
        self.header = self.data.columns
        
        self._log(
            'Finished compiling the `intercell network by resource` table.'
        )
    
    
    def setup_data(self):
        
        self.intercell = omnipath2.data.get_db('intercell')
        self.annot = omnipath2.data.get_db('annotations')
        self.network = omnipath2.data.get_db(self.network_dataset)
        self.network_df = omnipath2.data.network_df(self.network_dataset)
        self.intercell.register_network(self.network_df)
        
        self.intercell.set_interclass_network_df(
            only_class_levels = self.class_levels,
            only_proteins = self.only_proteins
        )
    
    
    def count_connections_pairwise(self):
        
        self._log('Counting connections between classes.')
        
        for mode in ('undirected', 'directed', 'stimulatory', 'inhibitory'):
            
            setattr(
                self,
                'con_%s' % mode,
                self.int_default(
                    getattr(
                        self.intercell,
                        'class_to_class_connections_%s' % mode
                    )()
                )
            )
        
        self.con_all = self.int_default(
            common.sum_dicts(
                self.con_undirected,
                self.con_directed,
                dict(
                    (
                        (cls1, cls0),
                        val
                    )
                    for (cls0, cls1), val in iteritems(self.con_directed)
                    if cls0 != cls1
                )
            )
        )
    
    
    def count_connections_groupwise(self):
        
        self._log('Counting degrees by class.')
        
        for mode, degrees_of in (('out', 'source'), ('in', 'target')):
            
            for con_mode in (
                'undirected',
                'directed',
                'stimulatory',
                'inhibitory',
            ):
                
                if mode == 'in' and con_mode == 'undirected':
                    
                    continue
                
                setattr(
                    self,
                    'degree_%s_%s' % (con_mode, mode),
                    self.int_default(
                        getattr(
                            self.intercell,
                            'degree_inter_class_network_%s_2' % con_mode
                        )(degrees_of = degrees_of)
                    )
                )
        
        self.degree_total = self.int_default(
            common.sum_dicts(
                self.degree_undirected_out,
                self.degree_directed_out,
                self.degree_directed_in,
            )
        )
    
    
    def count_connections(self):
    
        self._log('Counting connections in the network.')
        
        self.con_network = (
            len(self.network.interactions_undirected()) +
            len(self.network.interactions_directed())
        )
        self.con_network_undir = len(self.network.interactions_undirected())
        self.con_network_dir = len(self.network.interactions_directed())
        self.con_network_stim = len(self.network.interactions_stimulatory())
        self.con_network_inh = len(self.network.interactions_inhibitory())
    
    
    @staticmethod
    def int_default(d):
        
        return collections.defaultdict(int, d)


class AnnotationsByEntity(omnipath2.table.TableBase):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'fname': 'annots_by_entity_tsv',
        }
        param.update(kwargs)
        
        omnipath2.table.TableBase.__init__(self, **param)
    
    
    def load(self):
        
        self.complexdb = omnipath2.data.get_db('complex')
        self.annotdb = omnipath2.data.get_db('annotations')
        self.intercell = omnipath2.data.get_db('intercell')
        
        all_proteins = set(dataio.all_uniprots(swissprot = True))
        all_complexes = set(self.complexdb.complexes.keys())
        annotations = list(self.intercell.classes.keys())
        
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
        
        self.data = []
        
        for cls, elements in iteritems(self.intercell.classes):
            
            for elem in elements:
                
                self.data.append(
                    AnnotationRecord(
                        entity_id = elem.__str__(),
                        is_complex = entity.Entity._is_complex(elem),
                        cls = cls,
                        parent = self.intercell.parents[cls],
                        resource_label = (
                            self.intercell.get_resource_label(cls)
                        ),
                        class_label = self.intercell.get_class_label(cls),
                        class_type = self.intercell.get_class_type(cls),
                    )
                )
        
        self.data = pd.DataFrame(
            self.data,
            columns = AnnotationRecord._fields
        )
        self.header = self.data.columns


class ResourcesByEntity(omnipath2.table.TableBase):
    
    
    def __init__(self, network_dataset = 'omnipath', **kwargs):
        
        self.network_dataset = network_dataset
        
        param = {
            'fname': 'resources_by_entity_tsv',
            'fname_param': (
                self.network_dataset,
            )
        }
        param.update(kwargs)
        
        omnipath2.table.TableBase.__init__(self, **param)
    
    
    def load(self):
        
        ResourceRecord = collections.namedtuple(
            'ResourceRecord',
            [
                'entity_id',
                'resource',
                'is_complex',
                'resource_class',
            ],
        )
        
        self.network = omnipath2.data.network_df(self.network_dataset)
        
        all_entities = set(
            itertools.chain(
                *self.network.records[['id_a', 'id_b']].values
            )
        )
        
        entities_by_resource = self.network.entities_by_resource()
        
        self.data = []
        
        for resource, entities in iteritems(entities_by_resource):
            
            for entity_id in entities:
                
                self.data.append(
                    ResourceRecord(
                        entity_id = entity_id,
                        resource = resource,
                        is_complex = entity.Entity._is_complex(entity_id),
                        resource_class = (
                            data_formats.categories[resource]
                                if resource in data_formats.categories else
                            'z'
                        ),
                    )
                )
        
        self.data = pd.DataFrame(
            self.data,
            columns = ResourceRecord._fields,
        )
        self.header = self.data.columns
    
    
class ComplexesByResource(omnipath2.table.TableBase):
    
    
    def __init__(self, **kwargs):
        
        param = {
            'fname': 'complexes_by_resource_tsv',
        }
        param.update(kwargs)
        
        omnipath2.table.TableBase.__init__(self, **param)
    
    
    def load(self):
        
        ComplexRecord = collections.namedtuple(
            'ComplexRecord',
            [
                'complex_id',
                'resource',
            ],
        )
        
        self.complexdb = omnipath2.data.get_db('complex')
        
        self.data = []
        
        for cplex_name, cplex in iteritems(self.complexdb.complexes):
            
            for resource in cplex.sources:
                
                self.data.append(
                    ComplexRecord(
                        complex_id = cplex.__str__(),
                        resource = resource,
                    )
                )
        
        self.data = pd.DataFrame(self.data, columns = ComplexRecord._fields)
        self.header = self.data.columns


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
        
        self.collect_classes()
        self.export_intercell_coverages()
        self.export_intercell_coverages_by_resource()
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
