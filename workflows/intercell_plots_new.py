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

import itertools

import numpy as np
import scipy.cluster.hierarchy
import matplotlib as mpl
import matplotlib.cm
import pattern.en

from pypath import session_mod
from pypath import common

import data_tools
from data_tools.iterables import subsets
from data_tools.spatial import equidist_polar
from data_tools.iterables import similarity
from data_tools.plots import cluster_hmap
from data_tools.plots import upset_wrap
from data_tools.plots import chordplot

import workflows
from workflows import plot


class IntercellPlots(session_mod.Logger):
    
    
    def __init__(
            self,
            network_dataset = 'omnipath',
            **kwargs,
        ):
        
        session_mod.Logger.__init__(self, name = 'op2.intercell_plots')
        
        self.network_dataset = network_dataset
        
        for attr, val in kwargs.items():
            
            setattr(self, attr, val)
    
    
    def main(self):
        
        self.load()
        self.plot_ligands_per_receptor()
        self.plot_receptors_per_ligand()
        self.plot_protein_counts_by_class()
    
    
    def plot_ligands_per_receptor(self):
        
        param_attr = 'param_ligands_per_receptor_histo'
        
        param = getattr(self, param_attr) if hasattr(self, param_attr) else {}
        
        default_param = {
            'class0': 'ligand',
            'class1': 'receptor',
            'network_dataset': self.network_dataset,
        }
        default_param.update(param)
        param = default_param
        
        self.lig_per_rec_degree_histo = InterClassDegreeHisto(**param)
    
    
    def plot_receptors_per_ligand(self):
        
        param_attr = 'param_receptors_per_ligand_histo'
        
        param = getattr(self, param_attr) if hasattr(self, param_attr) else {}
        
        default_param = {
            'class0': 'ligand',
            'class1': 'receptor',
            'network_dataset': self.network_dataset,
        }
        default_param.update(param)
        param = default_param
        
        self.lig_per_rec_degree_histo = InterClassDegreeHisto(**param)
    
    
    def plot_protein_counts_by_class(self):
        
        self.counts_by_class_protein = CountsByClass()
    
    
    def plot_entity_counts_by_class(self):
        
        self.counts_by_class_all = CountsByClass(
            entity_types = {'protein', 'complex'}
        )
    
    
    def plot_class_similarities(self):
        
        self.csim = ClassSimilarities()


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
        
        self.data = workflows.data
        self.intercell = self.data.get_db('intercell')
        self.data.ensure_dataset(self.network_dataset)
        self.intercell.register_network(
            self.data.network_df(
                self.network_dataset,
                by_source = True,
            )
        )
        self.degrees = self.intercell.degree_inter_class_network(
            self.class0,
            self.class1,
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
            entity_types = 'protein',
            class_types = 'main',
            xscale_log = False,
            **kwargs,
        ):
        
        self.entity_types = common.to_set(entity_types)
        self.class_types = common.to_set(class_types)
        self.xscale_log = xscale_log
        
        param = {
            'maketitle': True,
            'xlab': 'Number of %s' % (
                pattern.en.pluralize(entity_types)
                    if isinstance(entity_types, common.basestring) else
                'entities'
            ),
            'legend': False,
        }
        param.update(kwargs)
        
        plot.PlotBase.__init__(self, **param)
    
    
    def load_data(self):
        
        self.labels, self.values = zip(*(
            sorted(
                self.counts.items(),
                key = lambda it: it[1],
                reverse = True,
            )
        ))
    
    
    def make_plots(self):
        
        self.get_subplot()
        self.ax.set_yticks(range(len(self.counts)))
        self.ax.grid()
        self.ax.scatter(self.values, range(len(self.counts)))
        self.ax.set_axisbelow(True)
        self.ax.set_yticklabels(labels = self.labels)
        self.ax.set_ylim(-1, len(self.counts))
        
        if self.xscale_log:
            
            self.ax.set_xscale('log')
        
        self.post_subplot_hook()


class CountsByClass(CountsScatterBase):
    
    
    def __init__(
        self,
        entity_types = 'protein',
        class_types = 'main',
        **kwargs,
    ):
        
        self.entity_types = common.to_set(entity_types)
        self.class_types = common.to_set(class_types)
        
        param = {
            'title': 'Entities by inter-cellular communication role',
            'fname': 'counts_by_class_pdf',
            'ylab': 'Inter-cellular communication roles',
            'fname_param': (
                '-'.join(sorted(self.entity_types)),
                '-'.join(sorted(self.class_types)),
            ),
        }
        param.update(kwargs)
        
        CountsScatterBase.__init__(
            self,
            entity_types = entity_types,
            class_types = class_types,
            **param
        )
    
    
    def load_data(self):
        
        self.data = workflows.data
        self.intercell = self.data.get_db('intercell')
        countsdf = self.intercell.counts_by_class(
            entity_types = self.entity_types,
            class_types = self.class_types,
        )
        self.counts = dict(zip(countsdf.index, countsdf))
        CountsScatterBase.load_data(self)


class CountsByResource(CountsScatterBase):
    
    
    def __init__(
        self,
        entity_types = 'protein',
        class_types = 'main',
        **kwargs,
    ):
        
        self.entity_types = common.to_set(entity_types)
        self.class_types = common.to_set(class_types)
        
        param = {
            'title': (
                'Entities in inter-cellular communication\n'
                'by annotation resource'
            ),
            'fname': 'counts_by_resource_pdf',
            'ylab': 'Inter-cellular annotation resources',
            'fname_param': (
                '-'.join(sorted(self.entity_types)),
            ),
            'height': 5,
        }
        param.update(kwargs)
        
        CountsScatterBase.__init__(
            self,
            entity_types = entity_types,
            class_types = class_types,
            **param
        )
    
    
    def load_data(self):
        
        self.data = workflows.data
        self.intercell = self.data.get_db('intercell')
        
        self.counts = self.intercell.counts_by_resource(
            entity_types = self.entity_types
        )
        CountsScatterBase.load_data(self)


class ClassSimilarities(plot.PlotBase):
    """
    Following ``data_tools.plots.cluster_hmap``.
    """
    
    
    def __init__(
            self,
            class_types = 'main',
            entity_types = 'protein',
            link_param = None,
            dendro_param = None,
            dendrogram_lwd = .5,
            **kwargs,
        ):
        
        self.entity_types = common.to_set(entity_types)
        self.class_types = common.to_set(class_types)
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
        
        self.data = workflows.data
        self.data.ensure_dataset('intercell')
        
        self.sims = []
        
        main_classes = sorted([
            cls
            for cls in self.data.intercell.classes.keys()
            if (
                cls in self.data.intercell.class_types and
                self.data.intercell.class_types[cls] == 'main'
            )
        ])
        
        for class0, class1 in itertools.product(main_classes, repeat = 2):
            
            self.sims.append(
                similarity(
                    self.data.intercell.classes[class0],
                    self.data.intercell.classes[class1],
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
        
        self.labels = np.array([
            self.data.intercell.class_labels[cls] for cls in main_classes
        ])
    
    
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
