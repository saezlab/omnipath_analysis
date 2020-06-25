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

import numpy as np

from pypath.share import session as session_mod


import omnipath2
import omnipath2.plot as plot


class CountsBase(plot.PlotBase):


    lightness_steps = (1.0, 1.33, 1.66, 2.0, 2.33)


    def __init__(
            self,
            network_dataset = 'omnipath',
            entity_type = 'protein',
            variables = (
                'entities',
                'interactions_0',
            ),
            variables2 = (
                (
                    'n_collection',
                    'n_shared_within_data_model',
                ),
                (
                    'n_collection',
                    'n_shared_within_data_model',
                ),
            ),
            by_category = True,
            split_by_categories = True,
            order_by = ('entities', 'n_collection'),
            palette = None,
            share_xaxis = True,
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
        self.share_xaxis = share_xaxis

        cols = len(variables2)

        param = {
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

        def get_value(var1, var2, key):

            collection = getattr(getattr(self, var1), var2)

            return collection[key] if key in collection else 0


        self.data = omnipath2.data
        self.network = self.data.get_db(self.network_dataset)

        for var in self.variables:

            setattr(
                self,
                var,
                getattr(
                    self.network,
                    'collect_%s' % var,
                )(
                    entity_type = self.entity_type,
                ),
            )

            this_var = getattr(self, var)

            del this_var.n_collection[('all', 'all', 'Total')]

            # adding totals by interaction type
            for s_u in ('shared', 'unique'):

                for k, v in getattr(
                    this_var,
                    'n_%s_by_interaction_type' % s_u
                ).items():

                    if k == 'all':

                        continue

                    getattr(
                        this_var,
                        'n_%s_within_data_model' % s_u
                    )[(k, 'all', 'Total')] = v

                for k, v in getattr(
                    this_var,
                    'n_%s_by_data_model' % s_u
                ).items():

                    if k[1] == 'all':

                        continue

                    getattr(
                        this_var,
                        'n_%s_within_interaction_type' % s_u
                    )[k + ('Total',)] = v


        self.main_var = getattr(self, self.order_by[0])

        self.keys = [
            it[0] for it in
            sorted(
                getattr(
                    self.main_var,
                    self.order_by[1]
                ).items(),
                key = lambda it: (
                    (
                        it[0][:2] # interaction type and data model
                            if self.by_category else
                        ()
                    ) +
                    (it[1],)
                        if it[0][2] != 'Total' else
                    ('ZZZ', 'ZZZ', 9999999)
                )
            )
        ]


        self.labels = np.array([
            (
                key[2]
                    if key[2] != 'Total' else
                key[1].replace('_', ' ').capitalize()
                    if key[1] != 'all' else
                key[0].replace('_', ' ').capitalize()
            )
            for key in self.keys
        ])

        if self.split_by_categories:

            self.cat_split = np.array([
                '__'.join(key[:2])
                for key in self.keys
            ])

        self.tick_loc = range(len(self.keys))

        self.values = [
            [
                np.array([
                    get_value(var1, var2, key)
                    for key in self.keys
                ])
                for var2 in self.variables2[ivar]
            ]
            for ivar, var1 in enumerate(self.variables)
        ]

        self.colors = [
            np.array([
                omnipath2.colors.lightness(
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
                key = lambda lab: 'ZZZ' if lab.endswith('all') else lab
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
                        # for scatterplot:
                        # facecolors = 'none',
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

                    self.ax.set_ylabel(
                        cat.split('__')[1].replace('_', ' ').capitalize()
                    )

                else:

                    self.ax.get_yaxis().label.set_visible(False)

                if self.share_xaxis and icat != self.grid_rows - 1:

                    self.ax.set_xticklabels([])
                    self.ax.get_xaxis().label.set_visible(False)

                else:

                    self.ax.set_xlabel(
                        'Number of %s' % (
                            ' '.join(
                                reversed(var1.split('_'))
                            ).replace('0 ', '')
                        )
                    )

        if self.share_xaxis:

            for ivar1 in range(len(self.variables)):

                xmins, xmaxs = zip(*(
                    ax[ivar1].get_xlim()
                    for ax in self.axes
                ))

                xmin = min(xmins)
                xmax = max(xmaxs)

                _ = list(
                    ax[ivar1].set_xlim(xmin, xmax)
                    for ax in self.axes
                )


class EdgeNodeCounts(CountsBase):


    def __init__(
        self,
        network_dataset = 'omnipath',
        share_xaxis = True,
        **kwargs
    ):

        param = {
            'variables': (
                'entities',
                'interactions_0',
            ),
            'variables2': (
                (
                    'n_collection',
                    'n_shared_within_interaction_type',
                    'n_shared_within_data_model',
                ),
                (
                    'n_collection',
                    'n_shared_within_interaction_type',
                    'n_shared_within_data_model',
                ),
            ),
            'fname': 'netw_node_edge_counts_pdf',
            'fname_param': (
                network_dataset,
                'shared-x'
                    if share_xaxis else
                'independent-x',
            ),
            'network_dataset': network_dataset,
            'share_xaxis': share_xaxis,
        }
        param.update(kwargs)

        CountsBase.__init__(self, **param)
