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


class InterClassDegreeHisto(plot.PlotBase):
    
    
    def __init__(
            self,
            class0,
            class1,
            label0 = None,
            label1 = None,
            network_dataset = 'omnipath',
            nbins = 100,
            **kwargs,
        ):
        
        self.network_dataset = network_dataset
        
        self.class0 = class0
        self.class1 = class1
        label0 = label0 or '%ss per %s' % (class0.capitalize(), class1)
        label1 = label1 or 'Frequency'
        self.nbins = nbins
        
        param = {
            'fname': 'inter_class_degree_pdf',
            'fname_param': (class0, class1),
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
                by_source = True
            )
        )
        self.degrees = self.intercell.degree_inter_class_network(
            self.class0,
            self.class1,
        )
    
    
    def make_plots(self):
        
        self.get_subplot()
        _ = self.plot_args.pop('cmap', None)
        self.ax.hist(self.degrees, bins = self.nbins, **self.plot_args)
        self.post_subplot_hook()
