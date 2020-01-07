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

import os
import subprocess

from pypath import session_mod

try:
    
    import PyPDF2 as pdflib
    
except ImportError:
    
    try:
        import pdfrw as pdflib
        
    except ImportError:
        
        session_mod.get_log().msg(
            'No PDF reader available (install either `PyPDF2` or `pdfrw`).',
            label = 'op2.latex',
        )

import omnipath2
import omnipath2.settings as op2_settings


class Latex(session_mod.Logger):
    
    
    def __init__(self):
        
        session_mod.Logger.__init__(self, name = 'op2.latex')
        
        self.latexdir = op2_settings.get('latex_dir')
        
        self.outdir = os.path.join(
            self.latexdir,
            op2_settings.get('latex_auto_dir'),
        )
    
    
    def main(self):
        
        self.collect_pdfs()
        self.write_inc_tex()
        self.write_all_fig_tex()
        self.latex_compile()
    
    
    def collect_pdfs(self):
        
        self._log('Collecting PDFs')
        
        self.pdfs = dict(
            (
                label,
                path
            )
            for label, path in
            omnipath2.files.files['recent'].items()
            if os.path.splitext(path) == '.pdf'
        )
        
        self._log('%s PDFs available.' % len(self.pdfs))
    
    
    def write_inc_tex(self):
        
        for label, pdfpath in self.pdfs.items():
            
            width, height = (
                pdflib.PdfFileReader(pdfpath).pages[0]['/MediaBox']
            )[-2:]
            
            size = (
                r'width=\textwidth'
                    if width > height * 1.1 else
                r'height=.8\textheight'
            )
            
            includegraphics = (
                r'\begin{figure}' '\n'
                '\t' r'\includegraphics[%s]{%s}\newpage' '\n'
                '\t' r'' '\n'
                r'\end{figure}' '\n'
                '\n' % (
                    size,
                    pdfpath,
                )
            )
            
            out_path = os.path.join(
                self.outdir,
                '%s.tex' % label,
            )
            
            with open(out_path, 'w') as fp:
                
                fp.write(includegraphics)
