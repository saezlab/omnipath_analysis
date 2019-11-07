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

import importlib as imp
import collections
import itertools

from pypath import session_mod

import omnipath2
from omnipath2 import settings as op2_settings
from omnipath2 import figures_preprocess
from omnipath2 import annotation_plots
from omnipath2 import intercell_plots
from omnipath2 import network_plots


_logger = session_mod.Logger(name = 'dataio')
_log = _logger._log


class ProductParam(object):
    
    
    def __init__(self, **kwargs):
        
        self.args, self.values = tuple(zip(*kwargs.items()))
    
    
    def __iter__(self):
        
        for val in itertools.product(*self.values):
            
            yield dict(zip(self.args, val))


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
    
    
    def run(self):
        
        _log('Running task `%s`.' % self.name)
        
        self.method(**self.param)
        
        _log('Task `%s` finished.' % self.name)


workflow = (
    
)


class Main(session_mod.Logger):
    
    
    def __init__(self):
        
        session_mod.Logger.__init__(self, name = 'op2.main')
        
        self._log(
            'OmniPath2 analysis and figures '
            'integrated workflow has been initialized.'
        )
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
        
        for obj in ('data', 'figpreproc'):
            
            if hasattr(self, obj):
                
                getattr(self, obj).reload()
    
    
    def main(self):
        
        self.load_data()
        self.data_for_r_plotting()
        self.run_r()
        self.make_annotation_plots()
        self.make_intercell_plots()
        self.make_network_plots()
    
    
    def load_data(self):
        
        self._log('Loading all databases.')
        
        self.data = omnipath2.data
        self.data.build()
    
    
    def data_for_r_plotting(self):
        
        self._log(
            'Creating `FiguresPreprocess` object '
            'for exporting tables for R plotting.'
        )
        self.figpreproc = figures_preprocess.FiguresPreprocess()
    
    
    def run_r(self):
        
        self._log('Running R plotting methods.')
        
        pass
    
    
    def make_annotation_plots(self):
        
        self.annotation_plots = annotation_plots_new.AnnotationPlots()
        self.annotation_plots.main()
    
    
    def make_intercell_plots(self):
        
        self.intercell_plots = intercell_plots_new.IntercellPlots()
        self.intercell_plots.main()
    
    
    def make_network_plots(self):
        
        self.network_plots = network_plots_new.NetworkPlots()
