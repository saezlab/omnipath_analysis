#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019-2020 Saez Lab
#
# OmniPath2 analysis and figures suite
#
# Authors:
#
# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de
#
# Dénes Türei
# turei.denes@gmail.com
#
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://omnipathdb.org/
#

from future.utils import iteritems

import importlib as imp
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

from pypath.inputs import main as dataio
from pypath.share import progress
from pypath.share import common
from pypath.core import entity

import omnipath2
from omnipath2 import settings as op2_settings
from omnipath2 import table


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
            'class_args': {
                'scope': 'generic',
                'source': 'composite',
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
                self.intercell.iter_classes(**self.class_args),
                2,
            )
                if self.mode == 'undirected' else
            itertools.product(
                self.intercell.iter_classes(**self.class_args),
                self.intercell.iter_classes(**self.class_args),
            )
        )

        for c0, c1 in iterator:

            numof_connections = getattr(
                self.intercell,
                method,
            )(
                annot_args_source = c0.args,
                annot_args_target = c1.args,
            )

            self.data.append([
                c0.name,
                c1.name,
                c0.name_label,
                c0.name_label,
                c0.n_proteins,
                c1.n_proteins,
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
                'cat_label',
                'category',
                'parent',
                'aspect',
                'scope',
                'source',
                'transmitter',
                'receiver',
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
            'protein': self.network.get_protein_identifiers(),
            'complex': self.network.get_complex_identifiers(),
        }

        self.data = []

        for cls in self.intercell.classes.values():

            for entity_type in ('protein', 'complex'):

                members = cls.filter_entity_type(entity_type = entity_type)
                total = len(members)
                in_network = len(
                    members &
                    network_entities[entity_type]
                )

                self.data.append([
                    cls.name_label,
                    cls.name,
                    cls.parent,
                    cls.aspect,
                    cls.scope,
                    cls.source,
                    cls.transmitter,
                    cls.receiver,
                    cls.resource,
                    entity_type,
                    total,
                    in_network,
                ])


class IntercellNetworkCounts(omnipath2.table.TableBase):


    def __init__(
            self,
            network_dataset = 'omnipath',
            only_proteins = True,
            annot_args = None,
            **kwargs
        ):

        self.network_dataset = network_dataset
        self.only_proteins = only_proteins
        self.annot_args = annot_args or {}
        self.entity_type = 'protein'

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

            a_cls0 = self.intercell.select(cls0, entity_type = entity_type)
            a_cls1 = self.intercell.select(cls1, entity_type = entity_type)

            a_parent0 = self.intercell.select(
                a_cls0.parent,
                entity_type = entity_type,
            )
            a_parent1 = self.intercell.select(
                a_cls1.parent,
                entity_type = entity_type,
            )

            overlap = a_cls0 & a_cls1

            in_network0 = in_network & a_cls0
            in_network1 = in_network & a_cls1

            cls01 = (a_cls0.name, a_cls1.name)
            cls10 = (a_cls1.name, a_cls0.name)

            self.data.append(
                ClassesPairwiseRecord(
                    name0 = a_cls0.name,
                    name1 = a_cls1.name,
                    label0 = a_cls0.name_label,
                    label1 = a_cls1.name_label,
                    resource0 = a_cls0.resource,
                    resource1 = a_cls1.resource,
                    aspect0 = a_cls0.aspect,
                    aspect1 = a_cls1.aspect,
                    source0 = a_cls0.source,
                    source1 = a_cls1.source,
                    scope0 = a_cls0.scope,
                    scope1 = a_cls1.scope,
                    transmitter0 = a_cls0.transmitter,
                    transmitter1 = a_cls1.transmitter,
                    receiver0 = a_cls0.receiver,
                    receiver1 = a_cls1.receiver,
                    total = len(reference_set), # total number of all proteins
                                                # or all complexes
                    network = len(in_network),
                    parent0 = a_parent0.name,
                    parent1 = a_parent1.name,
                    # sizes
                    size_cls0 = len(a_cls0),
                    size_cls1 = len(a_cls1),
                    overlap_cls01 = len(overlap),
                    in_network_cls0 = len(in_network0),
                    in_network_cls1 = len(in_network1),
                    size_parent0 = len(a_parent0),
                    size_parent1 = len(a_parent1),
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
                'name0',
                'name1',
                'label0',
                'label1',
                'resource0',
                'resource1',
                'aspect0',
                'aspect1',
                'source0',
                'source1',
                'scope0',
                'scope1',
                'transmitter0',
                'transmitter1',
                'receiver0',
                'receiver1',
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
        gc.collect()

        self.data = []

        annot_entities = (
            self.annot.proteins
                if self.only_proteins else
            self.annot.reference_set
        )

        self._log('Building the intercell network by resource table.')

        prg = progress.Progress(
            len(classes.items()) ** 2 / 2,
            'Collecting intercell categories stats',
            1,
            off = False,
        )

        entity_type = 'protein' if self.only_proteins else None
        in_network = self.network.get_identifiers(
            entity_type = self.entity_type
        )

        for cls0, cls1 in (
            itertools.combinations_with_replacement(
                self.intercell.filter_classes(
                    scope = 'generic',
                    aspect = 'functional',
                ),
                2
            )
        ):

            prg.step()

            add_stats_record(
                cls0.key,
                cls1.key,
                annot_entities,
            )

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
        """
        Ensures all required databases are loaded.
        """

        self.intercell = omnipath2.data.get_db('intercell')
        self.annot = omnipath2.data.get_db('annotations')
        self.network = omnipath2.data.get_db(self.network_dataset)
        self.network_df = omnipath2.data.network_df(self.network_dataset)
        self.intercell.register_network(self.network_df)


    def count_connections_pairwise(self):
        """
        Counts the network connections pairwise between the intercell classes.
        """

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
        """
        Counts the degrees for each of the intercell classes.
        """

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
        """
        Counts the connections in the entire network (not only the intercell
        annotated parts).
        """

        self._log('Counting connections in the network.')

        self.con_network = self.network.count_interactions_0()
        self.con_network_undir = (
            self.network.count_interactions_non_directed_0()
        )
        self.con_network_dir = self.network.count_interactions_directed()
        self.con_network_stim = self.network.count_interactions_positive()
        self.con_network_inh = self.network.count_interactions_negative()


    @staticmethod
    def int_default(d):

        return collections.defaultdict(int, d)


class IntercellAnnotationsByEntity(omnipath2.table.TableBase):


    def __init__(self, **kwargs):

        param = {
            'fname': 'intercell_annots_by_entity_tsv',
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
                'name',
                'label',
                'parent',
                'aspect',
                'source',
                'scope',
                'transmitter',
                'receiver',
                'resource',
            ],
        )

        self.data = []

        for key, cls in iteritems(self.intercell.classes):

            for elem in cls:

                self.data.append(
                    AnnotationRecord(
                        entity_id = elem.__str__(),
                        is_complex = entity.Entity._is_complex(elem),
                        name = cls.name,
                        labem = cls.name_label,
                        parent = cls.parent,
                        aspect = cls.aspect,
                        source = cls.source,
                        scope = cls.scope,
                        transmitter = cls.transmitter,
                        receiver = cls.receiver,
                        resource = cls.resource,
                    )
                )

        self.data = pd.DataFrame(
            self.data,
            columns = AnnotationRecord._fields
        )
        self.header = self.data.columns


class AnnotationsByEntity(omnipath2.table.TableBase):


    def __init__(self, **kwargs):

        param = {
            'fname': 'annots_by_entity_tsv',
        }
        param.update(kwargs)

        omnipath2.table.TableBase.__init__(self, **param)


    def load(self):

        self.annotdb = omnipath2.data.get_db('annotations')
        self.annotdb.make_narrow_df()

        self.data = self.annotdb.narrow_df


class EnzymeSubstrate(omnipath2.table.TableBase):


    def __init__(self, **kwargs):

        param = {
            'fname': 'enz_sub_tsv',
        }
        param.update(kwargs)

        omnipath2.table.TableBase.__init__(self, **param)


    def load(self):

        self.enzsubdb = omnipath2.data.get_db('enz_sub')
        self.enzsubdb.make_df()

        self.data = self.enzsubdb.df


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
                'interaction_type',
                'data_model',
            ],
        )

        self.network = omnipath2.data.get_db(self.network_dataset)

        all_entities = self.network.get_entities()

        entities_by_resource = (
            self.network.\
                entities_by_interaction_type_and_data_model_and_resource()
        )

        self.data = []

        for (itype, dmodel, res), entities in iteritems(entities_by_resource):

            for en in entities:

                self.data.append(
                    ResourceRecord(
                        entity_id = en.identifier,
                        resource = res,
                        is_complex = en.is_complex(),
                        interaction_type = itype,
                        data_model = dmodel,
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
                'references',
                'stoichiometry',
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
                        references = ';'.join(cplex.references),
                        stoichiometry = cplex.stoichiometry,
                    )
                )

        self.data = pd.DataFrame(
            self.data,
            columns = ComplexRecord._fields,
        )
        self.header = self.data.columns


class InterClassOverlaps(omnipath2.table.TableBase):


    def __init__(self, **kwargs):

        param = {
            'fname': 'category_overlaps_tsv',
            'header': [
                'cat0', 'cat1', 'size0', 'size1', 'total', 'overlap',
            ],
        }
        param.update(kwargs)

        omnipath2.table.TableBase.__init__(self, **param)


    def load(self):

        self.intercell = omnipath2.data.get_db('intercell')

        self.data = []

        for c0, c1 in (
            itertools.product(
                self.intercell.filter_classes(
                    scope = 'generic',
                    source = 'composite'
                ),
                self.intercell.filter_classes(
                    scope = 'generic',
                    source = 'composite'
                ),
            )
        ):

            self.data.append([
                c0.name,
                c1.name,
                len(c0.n_proteins),
                len(c1.n_proteins),
                len(c0.proteins | c1.proteins),
                len(c0.proteins & c1.proteins),
            ])


class NetworkCoverage(omnipath2.table.TableBase):


    groups = (
        ('kinase.com', 'Kinases'),
        ('TFcensus', 'TFs'),
        ('CancerGeneCensus', 'Cancer Drivers (COSMIC)'),
        ('IntOGen', 'Cancer Drivers (IntOGen)'),
        ('Phosphatome', 'Phosphatases'),
        ('HPMR', 'Receptors'),
    )


    def __init__(self, network_dataset = 'omnipath', **kwargs):

        param = {
            'fname': 'network_coverage_tsv',
            'fname_param': (network_dataset,),
            'header': [
                'resource',
                'resource_type',
                'group',
                'n_network',
                'n_group',
                'coverage',
            ],
        }
        param.update(kwargs)

        self.network_dataset = network_dataset

        omnipath2.table.TableBase.__init__(self, **param)


    def load(self):

        self.annot = omnipath2.data.get_db('annotations')
        self.network = omnipath2.data.get_db(self.network_dataset)

        self.data = []

        resources = self.network.get_resource_names()

        proteins = dict(
            (
                (resource, 'resource'),
                self.network.get_protein_identifiers(resources = resource)
            )
            for resource in resources
        )
        proteins[('OmniPath', 'total')] = (
            self.network.get_protein_identifiers()
        )
        proteins.update(
            dict(
                (
                    (data_model, 'data_model'),
                    self.network.get_protein_identifiers(
                        data_model = data_model
                    )
                )
                for data_model in self.network.get_data_models()
            )
        )

        for annot, label in self.groups:

            this_group = set(
                e
                for e in self.annot.annots[annot].annot.keys()
                if entity.Entity._is_protein(e)
            )

            for res, in_network in proteins.items():

                self.data.append([
                    (
                        res[0].capitalize().replace('_', ' ')
                            if res[1] == 'data_model' else
                        res[0]
                    ),
                    res[1],
                    label,
                    len(in_network),
                    len(this_group),
                    len(in_network & this_group),
                ])
