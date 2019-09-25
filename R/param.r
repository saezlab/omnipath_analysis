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

require(ggplot2)


get_default_param <- function(){
    
    list(
        
        typeface = 'HelveticaNeueLT Std Lt Cn',
        font_style = 'bold',
        dir_figures = file.path('figures', '20190925'),
        dir_data = file.path('tables', '20190925'),
        dir_misc = 'misc',
        add_titles = TRUE,
        theme = theme_minimal,
        theme_defaults = list(
            axis.text.x = element_text(color = '#000000'),
            axis.text.y = element_text(color = '#000000'),
            panel.grid.major.x = element_blank()
        ),
        input_intercell_cat_pairwise = 'stats_by_resource_20190925.tsv',
        input_intercell_annotation_by_entity = (
            'annotations_by_entity_20190814.tsv'
        ),
        input_resource_by_entity = 'resources_by_entity_20190925.tsv',
        input_complexes_by_resource = 'complexes_by_resource_20190925.tsv',
        graph_plot_defaults = list(
            vertex.frame.color = NA,
            vertex.color = '#FDC70F',
            edge.color = '#4C4B6B33',
            vertex.label.color = '#4C4B6B'
        ),
        graph_layout_default = 'fr',
        console_settings = list(
            width = 270,
            dplyr.width = 270
        ),
        width = 4,
        height = 3
    )
    
}
