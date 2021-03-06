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


import importlib as imp
import collections
import itertools

from pypath.share import common
from pypath.share import session as session_mod

import omnipath2
from omnipath2 import settings as op2_settings
from omnipath2 import r_preprocess
from omnipath2 import annotation_plots
from omnipath2 import intercell_plots
from omnipath2 import network_plots
from omnipath2 import complexes_plots
from omnipath2 import supptables
from omnipath2 import r_runner


_logger = session_mod.Logger(name = 'op2.main')
_log = _logger._log


class ProductParam(object):


    def __init__(self, **kwargs):

        self.args, self.values = tuple(zip(*kwargs.items()))


    def __iter__(self):

        for val in itertools.product(*self.values):

            yield dict(zip(self.args, val))


class IterParam(object):


    def __init__(self, *args):

        self.args = args


    def __iter__(self):

        for arg in self.args:

            yield arg


class Param(object):


    def __init__(self, **kwargs):

        self.param = kwargs


    def __iter__(self):

        yield self.param


class Task(
    collections.namedtuple(
        'TaskBase',
        [
            'method',
            'param',
            'name',
        ],
    )
):


    def __new__(cls, method, param = None, name = 'unknown'):

        param = param or Param()

        if not isinstance(param, (tuple, list)):

            param = param,

        return super(Task, cls).__new__(cls, method, param, name)


    def run(self):

        _log('Running task `%s`.' % self.name)

        for param in itertools.product(*self.param):

            param = dict(itertools.chain(*(par.items() for par in param)))

            _log(
                'Running with param `%s`.' % (
                    ', '.join(
                        '%s=%s' % (name, str(value))
                        for name, value in param.items()
                    )
                )
            )

            self.method(**param)

        _log('Task `%s` finished.' % self.name)


workflow = collections.OrderedDict(

    supptables = (
        Task(
            method = supptables.NetworkS2_PPIall,
            name = 'Supp Table S2, network all PPI',
        ),
        Task(
            method = supptables.NetworkS2_PPIcurated,
            name = 'Supp Table S2, network curated PPI',
        ),
        Task(
            method = supptables.NetworkS2_TFtarget,
            name = 'Supp Table S2, TF-target network',
        ),
        Task(
            method = supptables.NetworkS2_miRNAmRNA,
            name = 'Supp Table S2, miRNA-mRNA network',
        ),
        Task(
            method = supptables.NetworkS2_TFmiRNA,
            name = 'Supp Table S2, TF-miRNA network',
        ),
        Task(
            method = supptables.EnzSubS3,
            name = 'Supp Table S3, enzyme-substrate',
        ),
        Task(
            method = supptables.ComplexesS4,
            name = 'Supp Table S4, complexes',
        ),
        Task(
            method = supptables.AnnotationsS5,
            name = 'Supp Table S5, annotations',
        ),
        Task(
            method = supptables.IntercellS6,
            name = 'Supp Table S6, intercell',
        ),
    ),

    r_preprocess = (
        Task(
            method = r_preprocess.InterClassConnections,
            param = ProductParam(
                network_dataset = (
                    'omnipath',
                    'curated',
                ),
                mode = (
                    'all',
                    'undirected',
                    'directed',
                    'stimulatory',
                    'inhibitory',
                ),
            ),
            name = 'Interclass connections table',
        ),
        Task(
            method = r_preprocess.IntercellClasses,
            name = 'Intercell classes table',
        ),
        Task(
            method = r_preprocess.IntercellCoverages,
            param = ProductParam(
                network_dataset = (
                    'omnipath',
                    'curated',
                    'tf_target',
                ),
            ),
            name = 'Intercell network coverage table',
        ),
        Task(
            method = r_preprocess.IntercellNetworkCounts,
            param = ProductParam(
                network_dataset = (
                    'omnipath',
                    'curated',
                ),
            ),
            name = 'Intercell network counts table',
        ),
        Task(
            method = r_preprocess.IntercellAnnotationsByEntity,
            name = 'Intercell annotations by entity table',
        ),
        Task(
            method = r_preprocess.ComplexesByResource,
            name = 'Complexes by resource table',
        ),
        Task(
            method = r_preprocess.InterClassOverlaps,
            name = 'Intercell class overlaps table',
        ),
        Task(
            method = r_preprocess.IntercellNetworkByResource,
            name = 'Intercell network by resource',
        ),
        Task(
            method = r_preprocess.ResourcesByEntity,
            param = ProductParam(
                network_dataset = (
                    'omnipath',
                    'curated',
                    'tf_target',
                ),
            ),
            name = 'Resources by entity table',
        ),
        Task(
            method = r_preprocess.AnnotationsByEntity,
            name = 'Annotations by entity table',
        ),
        Task(
            method = r_preprocess.EnzymeSubstrate,
            name = 'Enzyme-substrate interactions table',
        ),
        Task(
            method = r_preprocess.NetworkCoverage,
            name = 'Network coverage on groups of proteins',
            param = ProductParam(
                network_dataset = (
                    'tf_target',
                    'omnipath',
                    'curated',
                ),
            ),
        ),
    ),

    complexes_plots = (
        Task(
            method = complexes_plots.ComplexesByResource,
            name = 'Complexes by resource figure',
        ),
    ),

    annotation_plots = (
        Task(
            method = annotation_plots.EntitiesByResource,
            name = 'Entities by annotation resource plot',
        ),
        Task(
            method = annotation_plots.RecordsByResource,
            name = 'Records by annotation resource plot',
        ),
        Task(
            method = annotation_plots.AnnotationNetworkOverlap,
            param = ProductParam(
                network_dataset = (
                    'omnipath',
                    'curated',
                    'tf_target',
                ),
            ),
            name = 'Annotation network overlap plot',
        ),
    ),

    network_plots = (
        Task(
            method = network_plots.EdgeNodeCounts,
            param = ProductParam(
                network_dataset = (
                    'omnipath',
                    'curated',
                    'tf_target',
                ),
                share_xaxis = (
                    True,
                    False,
                ),
            ),
            name = 'Network node and edge counts plot',
        ),
    ),

    intercell_plots = (
        Task(
            method = intercell_plots.InterClassDegreeHisto,
            param = (
                Param(
                    class0 = 'ligand',
                    class1 = 'receptor',
                ),
                ProductParam(
                    network_dataset = (
                        'omnipath',
                        'curated',
                    ),
                ),
            ),
            name = 'Ligand-receptor degrees histogram',
        ),
        Task(
            method = intercell_plots.CountsByClass,
            param = ProductParam(
                entity_type = (
                    'protein',
                    {'protein', 'complex'},
                ),
            ),
            name = 'Counts by intercell class plot',
        ),
        Task(
            method = intercell_plots.CountsByResource,
            param = ProductParam(
                entity_type = (
                    'protein',
                    {'protein', 'complex'},
                ),
            ),
            name = 'Counts by intercell resource plot',
        ),
        Task(
            method = intercell_plots.ClassSimilarities,
            name = 'Intercell class similarities plot',
        ),
        Task(
            method = intercell_plots.InterClassChordplot,
            param = ProductParam(
                network_dataset = (
                    'omnipath',
                    'curated',
                ),
                only_directed = (
                    True,
                    False,
                ),
            ),
            name = 'Intercell class connections chordplot',
        ),

    ),

    r_plotting = (
        Task(
            method = r_runner.RRunner,
            name = 'R plotting',
        ),
    ),

)


class Main(session_mod.Logger):


    def __init__(
            self,
            parts = None,
            steps = None,
        ):

        session_mod.Logger.__init__(self, name = 'op2.main')

        self.parts = common.to_set(parts)
        self.steps = (
            steps
                if isinstance(steps, (dict, type(None))) else
            {'main': steps}
        )

        self._log(
            'The OmniPath2 analysis and figures '
            'integrated workflow has been initialized.'
        )


    def reload(self):
        """
        Reloads the object from the module level.
        """

        omnipath2.data.reload()
        omnipath2.colors.reload()

        modules = (
            op2_settings,
            r_preprocess,
            network_plots,
            annotation_plots,
            complexes_plots,
            intercell_plots,
            supptables,
        )

        for mod in modules:

            imp.reload(mod)

        for dataset in omnipath2.data.datasets:

            if hasattr(omnipath2.data, dataset):

                getattr(omnipath2.data, dataset).reload()

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def main(self):

        self._log('Beginning workflow.')

        missing_parts = self.parts - set(workflow.keys())

        if missing_parts:

            self._log(
                'Warning: part(s) not defined in '
                'the workflow: %s.' % ', '.join(missing_parts)
            )

        _workflow = self.steps or workflow

        for part_name, part_tasks in _workflow.items():

            if self.parts and part_name not in self.parts:

                continue

            self._log('Beginning workflow part `%s`.' % part_name)

            for task in part_tasks:

                task.run()

        self._log('Workflow finished.')
