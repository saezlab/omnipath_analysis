#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019 Saez Lab
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

import subprocess

import pypath.session_mod

from omnipath2 import settings as op2_settings

class RRunner(pypath.session_mod.Logger):
    
    
    def __init__(
            self,
            r_command = None,
            **kwargs
        ):
        
        pypath.session_mod.Logger.__init__(self, name = 'op2.r_runner')
        
        self.r_command = r_command or op2_settings.get('r_workflow_command')
        
        self.run()
    
    
    def run(self):
        
        self._log('Calling `%s`.' % self.r_command)
        
        proc = subprocess.run(self.r_command.split())
        
        self.r_return_code = proc.returncode
        
        self._log(
            'Command `%s` finished with value %u.' % (
                self.r_command,
                self.r_return_code
            )
        )
