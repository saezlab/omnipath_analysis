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

import imp

from pypath import session_mod


import workflows
import workflows.plot as plot


class NetworkPlots(session_mod.Logger):
    
    
    def __init__(
            self,
            network_dataset = 'omnipath',
            **kwargs
        ):
        
        session_mod.Logger.__init__(self, name = 'op2.netw_plots')
        self._log('Compiling network plots.')
        
        self.network_dataset = network_dataset
    
    
    def main(self):
        
        pass


class CountsScatterBase(plot.PlotBase):
    
    
    lightness_steps = (1.0, 1.33, 1.66, 2.0, 2.33)
    
    
    def __init__(
            self,
            network_dataset = 'omnipath',
            entity_type = 'protein',
            variables = (
                'entities',
                'interactions_undirected',
            ),
            variables2 = (
                (
                    'by_resource',
                    'shared',
                    'shared_cat',
                ),
                (
                    'by_resource',
                    'shared',
                    'shared_cat',
                ),
            ),
            by_category = True,
            order_by = ('entities', 'by_resource'),
            palette = None,
            **kwargs
        ):
        
        self.network_dataset = network_dataset
        self.entity_type = entity_type
        self.variables = variables
        self.variables2 = variables2
        self.by_category = by_category
        self.order_by = order_by
        self.palette = palette
        
        cols = len(variables2)
        
        param = {
            'fname': 'netw_node_edge_counts_pdf',
            'fname_param': (self.network_dataset,),
            'ylab': 'Resources',
            'height': 7,
            'width': 3 * cols,
            'grid_cols': cols,
            'legend': False,
        }
        param.update(kwargs)
        
        plot.PlotBase.set_palette(self)
        plot.PlotBase.__init__(self, **param)
    
    
    def load_data(self):
        
        self.data = workflows.data
        self.network = self.data.get_db(self.network_dataset)
        
        for var in self.variables:
            
            setattr(
                self,
                var,
                getattr(
                    self.network,
                    '%s_stats' % var,
                )(
                    entity_type = self.entity_type,
                ),
            )
        
        self.main_var = getattr(self, self.order_by[0])
        
        self.labels = [
            it[0] for it in
            sorted(
                getattr(
                    self.main_var.counts,
                    self.order_by[1]
                ).items(),
                key = lambda it: (
                    (
                        (
                            self.main_var.counts.resource_cat[it[0]]
                                if it[0] in self.main_var.counts.resource_cat
                                else
                            it[0]
                        ,)
                            if self.by_category else
                        ()
                    ) +
                    (it[1],)
                )
            )
        ]
        
        self.tick_loc = range(len(self.labels))
        
        self.values = [
            [
                [
                    getattr(
                        getattr(self, var1).counts,
                        var2,
                    )[lab]
                    for lab in self.labels
                ]
                for var2 in self.variables2[ivar]
            ]
            for ivar, var1 in enumerate(self.variables)
        ]
        
        self.colors = [
            [
                workflows.colors.lightness(
                    self.palette.colors[ivar],
                    self.lightness_steps[ivar2],
                )
                for ivar2 in range(len(self.variables2[ivar]))
            ]
            for ivar, var1 in enumerate(self.variables)
        ]
    
    
    def make_plots(self):
        
        for ivar1, var1 in enumerate(self.variables):
            
            self.get_subplot(0, ivar1)
            
            for ivar2 in reversed(range(len(self.variables2[ivar1]))):
                
                self.ax.scatter(
                    self.values[ivar1][ivar2],
                    self.tick_loc,
                    color = self.colors[ivar1][ivar2],
                    facecolors = 'none',
                )
            
            self.ax.set_yticks(self.tick_loc)
            self.ax.set_yticklabels(self.labels)
            self.ax.grid()
            #self.ax.legend(loc = 0)
            self.ax.set_ylim(-1, len(self.tick_loc))
            self.ax.set_axisbelow(True)
