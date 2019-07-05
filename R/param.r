#!/usr/bin/env Rscript

#
#  This file is part of the `omnipath2` R package
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


get_default_param <- function(){
    
    list(
        
        typeface = 'DINPro',
        dir_figures = 'figures',
        dir_data = 'data',
        dir_misc = 'misc',
        add_titles = TRUE,
        theme = theme_bw
        
    )
    
}
