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
import collections

import numpy as np

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
                    'shared_res_cat',
                ),
                (
                    'by_resource',
                    'shared',
                    'shared_res_cat',
                ),
            ),
            by_category = True,
            split_by_categories = True,
            order_by = ('entities', 'by_resource'),
            palette = None,
            **kwargs
        ):
        
        self.network_dataset = network_dataset
        self.entity_type = entity_type
        self.variables = variables
        self.variables2 = variables2
        self.by_category = by_category
        self.split_by_categories = split_by_categories
        self.order_by = order_by
        self.palette = palette
        
        cols = len(variables2)
        
        param = {
            'fname': 'netw_node_edge_counts_pdf',
            'fname_param': (self.network_dataset,),
            'ylab': 'Resources',
            'height': 9,
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
        
        self.labels = np.array([
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
                        if it[0] != 'Total' else
                    ('ZZZ', 99999)
                )
            )
        ])
        
        if self.split_by_categories:
            
            self.cat_split = np.array([
                'Total'
                    if lab == 'Total' else
                self.main_var.counts.resource_cat[lab]
                for lab in self.labels
            ])
        
        self.tick_loc = range(len(self.labels))
        
        self.values = [
            [
                np.array([
                    getattr(
                        getattr(self, var1).counts,
                        var2,
                    )[lab]
                    for lab in self.labels
                ])
                for var2 in self.variables2[ivar]
            ]
            for ivar, var1 in enumerate(self.variables)
        ]
        
        self.colors = [
            np.array([
                workflows.colors.lightness(
                    self.palette.colors[ivar],
                    self.lightness_steps[ivar2],
                )
                for ivar2 in range(len(self.variables2[ivar]))
            ])
            for ivar, var1 in enumerate(self.variables)
        ]
        
        if self.split_by_categories:
            
            self.grid_rows = len(set(self.cat_split))
            cat_counts = collections.Counter(self.cat_split)
            self.categories = sorted(
                cat_counts.keys(),
                key = lambda lab: 'ZZZ' if lab == 'Total' else lab
            )
            self.grid_hratios = [cat_counts[cat] for cat in self.categories]
    
    
    def make_plots(self):
        
        for icat, cat in enumerate(
            self.categories
                if self.split_by_categories else
            ('Resources',)
        ):
            
            for ivar1, var1 in enumerate(self.variables):
                
                self.get_subplot(icat, ivar1)
                
                for ivar2 in range(len(self.variables2[ivar1])):
                    
                    tick_loc = (
                        range((self.cat_split == cat).sum())
                            if self.split_by_categories else
                        self.tick_loc
                    )
                    
                    self.ax.barh(
                        tick_loc,
                        self.values[ivar1][ivar2][
                            self.cat_split == cat
                                if self.split_by_categories else
                            np.array([True] * len(tick_loc))
                        ],
                        color = self.colors[ivar1][ivar2],
                        #facecolors = 'none',
                    )
                
                self.ax.set_yticks(tick_loc)
                
                if ivar1 == 0:
                    
                    self.ax.set_yticklabels(
                        self.labels[
                            self.cat_split == cat
                                if self.split_by_categories else
                            np.array([True] * len(tick_loc))
                        ]
                    )
                    
                else:
                    
                    self.ax.set_yticklabels([])
                
                self.ax.tick_params(axis = 'both', length = 0)
                self.ax.locator_params(nbins = 6, axis = 'x')
                
                self.ax.grid()
                #self.ax.legend(loc = 0)
                self.ax.set_ylim(-1, len(tick_loc))
                self.ax.set_axisbelow(True)
                
                self.post_subplot_hook()
                
                if ivar1 == 0:
                    
                    self.ax.set_ylabel(cat)
                    
                else:
                    
                    self.ax.get_yaxis().label.set_visible(False)
                
                if icat != self.grid_rows - 1:
                    
                    self.ax.set_xticklabels([])
                    self.ax.get_xaxis().label.set_visible(False)
                    
                else:
                    
                    self.ax.set_xlabel(
                        'Number of %s' % var1.replace('_', ' ')
                    )
