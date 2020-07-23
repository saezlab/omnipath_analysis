#!/usr/bin/env Rscript

#
#  This file is part of the `omnipath2` R package
#
#  Copyright
#  2019-2020
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
            axis.text.y = element_text(color = '#000000')
            #panel.grid.minor.x = element_blank()
        ),
        input_intercell_cat_pairwise_tsv = 'stats_by_resource_%s_%s',
        input_intercell_annotation_by_entity_tsv = {
            'intercell_annotations_by_entity'
        },
        input_annotation_by_entity_tsv = 'annotations_by_entity',
        input_enzyme_substrate_tsv = 'enzyme_substrate',
        input_resource_by_entity_tsv = 'resources_by_entity_%s',
        input_complexes_by_resource_tsv = 'complexes_by_resource_tsv',
        input_network_coverages_tsv = 'network_coverage_%s',
        input_network_size_tsv = 'S2_network_%s_for-r',
        input_intercell_network_by_resource_tsv = {
            'intercell_network_by_resource'
        },
        input_intercell_network_summary_tsv = (
            'inter_class_network_summary_%s_%s_%s'
        ),
        
        graph_plot_defaults = list(
            vertex.frame.color = NA,
            vertex.color = '#FDC70F',
            edge.color = '#4C4B6B33',
            vertex.label.color = '#4C4B6B'
        ),
        graph_layout_default = 'fr',
        
        two_shades_1 = c('#4268B3', '#B3C5E9'),
        three_shades_1 = c('#4268B3', '#6F8DCF', '#B3C5E9'),
        palette = c('#176FC1', '#00AAB0'),
        palette2 = c(
            '#176FC1',
            '#A6D81C',
            '#F89D0E',
            '#ED0772',
            '#0A6167',
            '#00AAB0',
            '#D22027',
            '#9E1639',
            '#7264B9',
            '#5B205F'
        ),
        
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

        intercell_main_classes = list(
            'ecm',
            'adhesion',
            'ligand',
            'receptor',
            'cell_surface_enzyme',
            'secreted_enzyme',
            'transporter'
        ),

        network_datasets = list(
            'omnipath',
            'tf_target',
            'mirna_mrna',
            'tf_mirna'
        ),

        network_dataset_labels = list(
            omnipath = 'post-translational',
            tf_target = 'transcriptional',
            mirna_mrna = 'post-transcriptional',
            tf_mirna = 'miRNA transcriptional'
        ),

        intercell_classes_exclude = list(
            'Small molecule ligand synthase',
            'Ligand dgidb',
            'Receptor regulator',
            'Matrix adhesion regulator',
            'Intracellular intercellular related',
            'Ligand regulator',
            'Ecm regulator',
            'Desmosome',
            'Matrix adhesion regulator',
            'Ion channel regulator',
            'Secreted peptidase inhibitor',
            'Cell surface peptidase',
            'Secreted receptor'
        ),

        intercell_colors = list(
            `Adhesion` = '#7264B9',
            `Adherens junction` = '#7264B9',
            `Tight junction` = '#7264B9',
            `Matrix adhesion` = '#7264B9',
            `Cell adhesion` = '#7264B9',
            `Secreted enzyme` = '#9E1639',
            `Secreted peptidase` = '#9E1639',
            `Ecm` = '#591E5E',
            `Ligand` = '#F89D0E',
            `Receptor` = '#00AAB0',
            `Cell surface enzyme` = '#ED0772',
            `Cell surface ligand` = '#D22027',
            `Transporter` = '#FFF200',
            `Gap junction` = '#FFF200',
            `Ion channel` = '#FFF200'
        ),

        intercell_colors_2 = list(
            `Adhesion` = '#CDC4E7',
            `Adherens junction` = '#CDC4E7',
            `Tight junction` = '#CDC4E7',
            `Matrix adhesion` = '#CDC4E7',
            `Cell adhesion` = '#CDC4E7',
            `Secreted enzyme` = '#8B456E',
            `Secreted peptidase` = '#8B456E',
            `Ecm` = '#BEA2BE',
            `Ligand` = '#FDDCA6',
            `Receptor` = '#ABE3DF',
            `Cell surface enzyme` = '#F9ABCA',
            `Cell surface ligand` = '#E16C70',
            `Transporter` = '#FFFBAB',
            `Gap junction` = '#FFFBAB',
            `Ion channel` = '#FFFBAB'
        ),

        fig_res_cov = 'res-cov_%s_%s_%s',
        fig_con_enrich = 'connection-enrichment_%s_%s_%s_%s_%s',
        fig_subcls_intersect = 'subclass-intersection_%s_%s_%s',
        fig_cplex_overlap = 'complex-overlap',
        fig_intercell_cls_graph = 'intercell-classes-graph_%s_%s',
        fig_cat_sizes = 'category-sizes_%s_%s_%s_%s',
        fig_cat_sizes_dots = 'category-sizes-dots',
        fig_enzyme_substrate_shared = 'enzyme-substrate-shared',
        fig_enzyme_substrate_self = 'enzyme-substrate-self-cplex-other_%s_%s',
        fig_enzyme_substrate_modtype = 'enzyme-substrate-modtype',
        fig_enzyme_substrate_modtype_dot = 'enzyme-substrate-modtype-dot',
        fig_enzyme_substrate_numof_res = 'enzyme-substrate-numof-res',
        fig_enzyme_by_substrate = 'enzyme-by-substrate',
        fig_substrate_by_enzyme = 'substrate-by-enzyme',
        fig_substrate_numof_sites = 'substrate-numof-sites',
        fig_enzyme_substrate_numof_ref = 'enzyme-substrate-numof-ref',
        fig_enzyme_substrate_by_ref = 'enzyme-substrate-by-ref',
        fig_cplex_by_resource = 'complexes-by-resource_%s_%s',
        fig_cplex_n_comp = 'complexes-numof-components',
        fig_refs_by_cplex = 'references-by-complex',
        fig_complexes_by_ref = 'complexes-by-ref',
        fig_cplex_by_comp = 'complexes-by-components',
        fig_network_coverage = 'network-coverage_%s_%s',
        fig_network_size = 'network-size_%s_%s',
        fig_network_dir = 'network-dir_%s_%s',
        fig_intercell_deg = 'intercell-degrees_%s_%s_%s_%s',
        fig_intercell_cov = 'intercell-cov_%s_%s_%s_%s',
        fig_intercell_size = 'intercell-size_%s_%s_%s_%s',
        fig_intercell_by_res = 'intercell-by-resource',
        fig_interc_chordplot = 'intercell-chordplot_%s_%s_%s'
    )
    
}
