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


class NetworkPlots(session_mod.Logger):
    
    
    def __init__(self, network_dataset = 'omnipath'):
        
        session_mod.Logger.__init__(self, name = 'op2.netw_plots')
        self._log('Compiling network plots.')
        
        self.network_dataset = network_dataset
    
    
    def main(self):
        
        pass
