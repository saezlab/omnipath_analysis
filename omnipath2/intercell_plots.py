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

import itertools

import numpy as np
import scipy.cluster.hierarchy
import matplotlib as mpl
import matplotlib.cm
import pattern.en

from pypath.share import session as session_mod
from pypath.share import common

import data_tools
from data_tools.iterables import subsets
from data_tools.spatial import equidist_polar
from data_tools.iterables import similarity
from data_tools.plots import cluster_hmap
from data_tools.plots import upset_wrap
from data_tools.plots import chordplot

import omnipath2
from omnipath2 import plot
from omnipath2 import settings


class InterClassDegreeHisto(plot.PlotBase):


    def __init__(
            self,
            class0,
            class1,
            label0 = None,
            label1 = None,
            network_dataset = 'omnipath',
            only_directed = True,
            only_effect = None,
            degrees_of = 'target',
            nbins = 100,
            log_y = False,
            **kwargs,
        ):

        self.network_dataset = network_dataset

        self.class0 = class0
        self.class1 = class1
        _class0, _class1 = (
            (class0, class1)
                if degrees_of == 'target' else
            (class1, class0)
        )
        label0 = label0 or '%ss per %s' % (_class0.capitalize(), _class1)
        label1 = label1 or 'Frequency'
        self.nbins = nbins
        self.only_directed = only_directed
        self.only_effect = only_effect
        self.log_y = log_y
        self.degrees_of = degrees_of

        param = {
            'fname': 'inter_class_degree_pdf',
            'fname_param': (
                class0 if self.degrees_of == 'target' else class1,
                class1 if self.degrees_of == 'target' else class0,
                'directed' if only_directed else 'undirected',
                (
                    'stimulation'
                        if only_effect == 1 else
                    'inhibition'
                        if only_effect == -1 else
                    'any_effect'
                ),
            ),
            'xlab': label0,
            'ylab': label1,
            'title': label0,
            'legend': False,
        }
        param.update(kwargs)

        plot.PlotBase.__init__(self, **param)


    def load_data(self):

        self.data = omnipath2.data
        self.intercell = self.data.get_db('intercell')
        self.data.ensure_dataset(self.network_dataset)
        self.intercell.register_network(
            self.data.network_df(
                self.network_dataset,
                by_source = True,
            )
        )
        self.degrees = self.intercell.degree_inter_class_network(
            annot_args_source = {
                'name': self.class0,
                'source': 'composite',
            },
            annot_args_target = {
                'name': self.class1,
                'source': 'composite',
            },
            only_directed = self.only_directed,
            only_effect = self.only_effect,
            degrees_of = self.degrees_of,
        )


    def make_plots(self):

        self.get_subplot()
        _ = self.plot_args.pop('cmap', None)
        self.ax.hist(self.degrees, bins = self.nbins, **self.plot_args)

        if self.log_y:

            self.ax.set_yscale('log')

        self.post_subplot_hook()


class CountsScatterBase(plot.PlotBase):


    def __init__(
            self,
            entity_type = 'protein',
            annot_args = None,
            xscale_log = False,
            **kwargs,
        ):

        self.annot_args = annot_args or {
            'source': 'composite',
        }
        self.entity_type = common.to_set(entity_type)
        self.xscale_log = xscale_log

        param = {
            'maketitle': True,
            'xlab': 'Number of %s' % (
                pattern.en.pluralize(entity_type)
                    if isinstance(entity_type, common.basestring) else
                'entities'
            ),
            'legend': False,
        }
        param.update(kwargs)

        plot.PlotBase.__init__(self, **param)


    def load_data(self):

        self.counts_nonzero = dict(
            it
            for it in self.counts.items()
            if it[1] != 0
        )

        self.labels, self.values = zip(*(
            sorted(
                self.counts_nonzero.items(),
                key = lambda it: it[1],
                reverse = True,
            )
        ))


    def make_plots(self):

        self.get_subplot()
        self.ax.set_yticks(range(len(self.counts_nonzero)))
        self.ax.grid()
        self.ax.scatter(self.values, range(len(self.counts_nonzero)))
        self.ax.set_axisbelow(True)
        self.ax.set_yticklabels(labels = self.labels)
        self.ax.set_ylim(-1, len(self.counts_nonzero))

        if self.xscale_log:

            self.ax.set_xscale('log')

        self.post_subplot_hook()


class CountsByClass(CountsScatterBase):


    def __init__(
        self,
        entity_type = 'protein',
        annot_args = None,
        **kwargs,
    ):

        self.entity_type = common.to_set(entity_type)

        param = {
            'title': 'Entities by inter-cellular communication role',
            'fname': 'counts_by_class_pdf',
            'ylab': 'Inter-cellular communication roles',
            'fname_param': (
                '-'.join(sorted(self.entity_type)),
            ),
        }
        param.update(kwargs)

        CountsScatterBase.__init__(
            self,
            entity_type = entity_type,
            annot_args = annot_args,
            **param
        )


    def load_data(self):

        self.data = omnipath2.data
        self.intercell = self.data.get_db('intercell')
        countsdf = self.intercell.counts_df(
            entity_type = self.entity_type,
            **self.annot_args
        )
        self.counts = dict(zip(
            (cat.capitalize() for cat in countsdf.category),
            countsdf.n_uniprot
        ))
        CountsScatterBase.load_data(self)


class CountsByResource(CountsScatterBase):


    def __init__(
        self,
        entity_type = 'protein',
        annot_args = None,
        **kwargs,
    ):

        self.entity_type = common.to_set(entity_type)

        param = {
            'title': (
                'Entities in inter-cellular communication\n'
                'by annotation resource'
            ),
            'fname': 'counts_by_resource_pdf',
            'ylab': 'Inter-cellular annotation resources',
            'fname_param': (
                '-'.join(sorted(self.entity_type)),
            ),
            'height': 5,
        }
        param.update(kwargs)

        CountsScatterBase.__init__(
            self,
            entity_type = entity_type,
            annot_args = annot_args,
            **param
        )


    def load_data(self):

        self.data = omnipath2.data
        self.intercell = self.data.get_db('intercell')

        self.counts = self.intercell.counts_by_resource(
            entity_types = self.entity_type
        )
        CountsScatterBase.load_data(self)


class ClassSimilarities(plot.PlotBase):
    """
    Following ``data_tools.plots.cluster_hmap``.
    """


    def __init__(
            self,
            annot_args = None,
            entity_types = 'protein',
            link_param = None,
            dendro_param = None,
            dendrogram_lwd = .5,
            **kwargs,
        ):

        self.entity_types = common.to_set(entity_types)
        self.annot_args = annot_args or {
            'source': 'composite',
        }
        self.link_param = link_param or {}
        self.dendro_param = dendro_param or {}
        self.dendrogram_lwd = dendrogram_lwd

        param = {
            'maketitle': True,
            'grid_rows': 2,
            'grid_cols': 3,
            'grid_hratios': (1, 7),
            'grid_wratios': (7, 1, 1),
            'width': 4,
            'height': 4,
            'palette': mpl.cm.viridis,
            'fname': 'inter_class_sim_pdf',
            'title': (
                'Szymkiewicz–Simpson similarity\n'
                'of major intercellular classes'
            ),
            'ylab': 'Inter-cellular communication roles',
            'xlab': 'Inter-cellular communication roles',
            'tight_layout': True,
        }
        param.update(kwargs)

        plot.PlotBase.__init__(self, **param)


    def load_data(self):

        self.data = omnipath2.data
        self.data.ensure_dataset('intercell')

        self.sims = []

        main_classes = sorted(
            (
                cls
                for cls in self.data.intercell.iter_classes(
                    **self.annot_args
                )
            ),
            key = lambda c: c.name_label,
        )

        for class0, class1 in itertools.product(main_classes, repeat = 2):

            self.sims.append(
                similarity(
                    class0,
                    class1,
                    mode = 'ss'
                )
            )

        self.sims = np.array(self.sims).reshape(
            len(main_classes),
            len(main_classes),
        )
        self.xlinked = scipy.cluster.hierarchy.linkage(
            self.sims.T,
            **self.link_param
        )
        self.ylinked = scipy.cluster.hierarchy.linkage(
            self.sims,
            **self.link_param
        )

        self.labels = np.array([cls.name_label for cls in main_classes])


    def make_plots(self):

        self.make_dendrogram_right()
        self.make_dendrogram_top()
        self.make_heatmap()

        self.axes[1][0].get_shared_x_axes().join(
            self.axes[0][0],
            self.axes[1][0],
        )
        self.axes[1][0].get_shared_y_axes().join(
            self.axes[1][1],
            self.axes[1][0],
        )

        self.get_subplot(1, 2)
        self.colorbar = self.fig.colorbar(
            self.heatmap,
            cax = self.ax,
            ax = self.ax,
            orientation = 'vertical',
        )
        self.set_title()
        self.fig.subplots_adjust(hspace = 0)
        self.fig.subplots_adjust(wspace = 0)
        self.colorbar.outline.set_visible(False)



    def make_heatmap(self):

        sims_ordered = self.sims[
            :,
            self.top_dendro['leaves']
        ][
            self.right_dendro['leaves'],
            :
        ]
        xlabels_ordered = self.labels[self.top_dendro['leaves']]
        ylabels_ordered = self.labels[self.right_dendro['leaves']]


        self.get_subplot(1, 0)

        self.heatmap = self.ax.imshow(
            sims_ordered,
            interpolation = 'none',
            cmap = self.palette,
            aspect = 'auto',
        )

        self.ax.set_xticks(range(len(xlabels_ordered)))
        self.ax.set_yticks(range(len(ylabels_ordered)))

        self.ax.set_xticklabels(xlabels_ordered, rotation = 90)
        self.ax.set_yticklabels(ylabels_ordered)
        ylim = self.ax.get_ylim()
        self.ax.set_ylim(ylim[0] + .5, ylim[1] - .5)
        self.ax.tick_params(axis = 'both', which = 'both', length = 0)
        [axis.set_linewidth(0) for axis in self.ax.spines.values()]


    def make_dendrogram_top(self):

        self._make_dendrogram(0, 0)


    def make_dendrogram_right(self):

        self._make_dendrogram(1, 1, orientation = 'right')


    def _make_dendrogram(self, i, j, orientation = 'top'):

        dendro_param = {}
        dendro_param.update(self.dendro_param)

        self.get_subplot(i, j)

        with mpl.rc_context({'lines.linewidth': self.dendrogram_lwd}):

            setattr(
                self,
                '%s_dendro' % orientation,
                scipy.cluster.hierarchy.dendrogram(
                    self.xlinked if orientation == 'top' else self.ylinked,
                    ax = self.ax,
                    orientation = orientation,
                    link_color_func = lambda k: 'k',
                    **dendro_param
                )
            )

        if orientation == 'right':

            self.ax.set_ylim(10 * self.sims.shape[0], 0)

        self.ax.set_axis_off()


class InterClassChordplot(plot.PlotBase):


    def __init__(
            self,
            network_dataset = 'omnipath',
            intercell_network_param = None,
            annot_args = None,
            only_directed = True,
            **kwargs
        ):

        self.network_dataset = network_dataset
        self.intercell_main_classes = settings.get('intercell_main_classes')

        self.annot_args = {
            'source': 'composite',
            'aspect': 'functional',
            'category': self.intercell_main_classes
        }
        self.annot_args.update(annot_args or {})

        self.intercell_network_param = {
            'annot_args': self.annot_args,
            'only_directed': only_directed,
            'only_effect': None,
        }
        self.intercell_network_param.update(intercell_network_param or {})

        icnparam = self.intercell_network_param

        param = {
            'fname': 'inter_class_chordplot_pdf',
            'fname_param': (
                self.network_dataset,
                'directed'
                    if icnparam['only_directed'] else
                'undirected',
                (
                    'stimulation'
                        if icnparam['only_effect'] == 1 else
                    'inhibition'
                        if icnparam['only_effect'] == -1 else
                    'any_effect'
                ),
            ),
            'palette': omnipath2.colors.get_palette('rwth_selected.gpl'),
            'make_plot_first': True,
            'legend': False,
            'legend_font_size': 8,
            'width': 10,
            'height': 8,
        }
        param.update(kwargs)

        plot.PlotBase.__init__(self, **param)


    def load_data(self):

        self.data = omnipath2.data
        self.intercell = self.data.get_db('intercell')
        network = self.data.network_df(self.network_dataset)
        self.intercell.register_network(network)
        self.edges = self.intercell.class_to_class_connections(
            **self.intercell_network_param,
        )
        self.edges = self.edges.reindex(
            index = itertools.product(*((self.intercell_main_classes,) * 2))
        )
        self.edges.index.set_levels(
            [
                [
                    cls.capitalize().replace('_', ' ')
                    for cls in self.edges.index.levels[0]
                ],
                [
                    cls.capitalize().replace('_', ' ')
                    for cls in self.edges.index.levels[1]
                ],
            ],
            inplace = True,
        )

        self.edges.rename('connections')
        self.edges = (
            self.edges.reset_index().rename(columns = {0: 'connections'})
        )

        self.adjacency = self.edges.pivot(
            index = 'category_a',
            columns = 'category_b',
            values = 'connections',
        )
        counts = self.intercell.counts_df(**self.annot_args)
        counts['label'] = [
            c.capitalize().replace('_', ' ')
            for c in counts.category
        ]
        counts.set_index('label', inplace = True)
        self.segments = counts.n_uniprot
        self.segments.name = 'size'
        self.segments = self.segments.reindex(
            index = self.edges.category_a.unique()
        )
        self.labels = self.segments.index


    def make_plots(self):

        self.colors = [
            self.palette.colors[i]
            for i in range(len(self.segments))
        ]

        self.fig = data_tools.plots.chordplot(
            edges = self.edges,
            nodes = self.segments,
            colors = self.colors,
            labels = False,
            alpha = .5,
        )

        self.ax = self.fig.gca()
        self.ax.set_position([0, 0, 0.75, 1]) # [left, bottom, width, height]
                                              # chordplot

        ax2 = self.fig.add_axes([0.70, 0.75, 0.25, 0.25]) # legend

        ax2.legend(
            [
                matplotlib.lines.Line2D(
                    [0],
                    [0],
                    marker = 'o',
                    color = 'w',
                    markerfacecolor = c
                )
                for c in self.colors
            ],
            self.labels,
            loc = 'center',
            ncol = 2,
        )

        ax2.set_axis_off()

        # barplot
        ax3 = self.fig.add_axes([0.72, 0.5, 0.25, 0.25])
        # heatmap
        ax4 = self.fig.add_axes([0.72, .05, 0.25, 0.5], sharex = ax3)

        ax3.bar(
            range(len(self.segments)),
            self.segments.values,
            color = self.colors,
        )
        #ax3.set_axis_off()
        ax3.spines['right'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax3.spines['bottom'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax3.set_xticks([])
        #ax.axes.get_yaxis().set_visible(False)
        ax3.locator_params(axis = 'y', nbins = 4)

        #ax.tick_params(axis='both', left='off', top='off', right='off',
        # bottom='off', labelleft='off', labeltop='off', labelright='off',
        # labelbottom='off')

        #ax3.get_yaxis().set_visible(True)

        adj = np.triu(self.adjacency.values).astype(float)
        adj[np.where(adj == 0)] = np.nan

        self.heatmap = ax4.imshow(adj, cmap = 'plasma')
        # [left, bottom, width, height]
        ax5 = self.fig.add_axes([0.72, 0.05, 0.25, 0.03]) # colorbar
        colorbar = self.fig.colorbar(
            self.heatmap,
            cax = ax5,
            ax = ax5,
            orientation = 'horizontal',
        )
        colorbar.outline.set_visible(False)
        ax4.set_axis_off()


    def post_subplot_hook(self):

        pass

    def finish(self):

        self.save()