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


class CountsByClass(plot.PlotBase):
    
    
    def __init__(self, entity_types = 'protein', class_types = 'main'):
        
        self.entity_types = common.to_set(entity_types)
        self.class_types = common.to_set(class_types)
        
        plot.PlotBase.__init__(
            self,
            fname = 'counts_by_class_pdf',
            fname_param = (
                '-'.join(sorted(self.entity_types)),
                '-'.join(sorted(self.class_types)),
            ),
            title = 'Entities by inter-cellular communication role',
            xlab = 'Number of entitites',
            ylab = 'Inter-cellular communication roles',
            maketitle = True,
        )
    
    
    def load_data(self):
        
        self.data = workflows.data
        self.intercell = self.data.get_db('intercell')
        self.counts = self.intercell.counts_by_class(
            entity_types = self.entity_types,
            class_types = self.class_types,
        )
    
    
    def make_plots(self):
        
        self.get_subplot()
        self.ax.set_yticks(range(len(self.counts)))
        self.ax.scatter(self.counts, range(len(self.counts)))
        self.ax.set_yticklabels(labels = self.counts.index)
        self.ax.set_ylim(-1, len(self.counts))
        self.ax.grid()
        #self.ax.set_xscale('log')


class ClassSimilarities(plot.PlotBase):
    
    
    def __init__(self, class_types = 'main', entity_types = 'protein'):
        
        self.entity_types = common.to_set(entity_types)
        self.class_types = common.to_set(class_types)
        
        plot.PlotBase.__init__(self)
        
    
    def load_data(self):
        
        sims = []
        
        for (class0, class1) in itertools.product(
                self.data.intercell.classes.keys(),
                repeat = 2,
            ):
            
            sims.append(similarity(set(a), set(b), mode='ss'))

        sims = np.array(sims).reshape(len(groups), len(groups))
