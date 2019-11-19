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
        files_json = 'files.json',
        dir_figures = 'figures',
        dir_tables = 'tables',
        dir_misc = 'misc',
        dir_log = 'omnipath2_log',
        add_titles = TRUE,
        theme = theme_minimal,
        theme_defaults = list(
            axis.text.x = element_text(color = '#000000'),
            axis.text.y = element_text(color = '#000000'),
            panel.grid.major.x = element_blank()
        ),
        input_intercell_cat_pairwise_tsv = 'stats_by_resource_%s_%s_%s',
        input_intercell_annotation_by_entity_tsv = 'annotations_by_entity',
        input_resource_by_entity_tsv = 'resources_by_entity_%s',
        input_complexes_by_resource_tsv = 'complexes_by_resource',
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
        height = 3,
        timestamp_format = '%Y%m%d',
        timestamp_dirs = TRUE,
        timestamp_files = TRUE,
        loglevel = 5,
        
        fig_res_cov = 'res-cov_%s_%s_%s_%s',
        fig_con_enrich = 'connection-enrichment_%s_%s_%s_%s',
        fig_subcls_intersect = 'subclass-intersection_%s_%s_%s',
        fig_cplex_overlap = 'complex-overlap',
        fig_intercell_cls_graph = 'intercell-classes-graph_%s_%s_%s',
        fig_cat_sizes = 'category-sizes_%s_%s_%s'
        
    )
    
}
