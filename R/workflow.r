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


con_enrich_input_param <- ProductParam$new(
    list('curated', 'omnipath'),
    list('proteins')
)

res_by_entity_input_param <- Param$new('omnipath', 'curated')

network_param <- Param$new(
    'omnipath',
    'curated',
    'tf_target',
    'mirna_mrna',
    'tf_mirna'
)


omnipath2_workflow <<- list(
    
    # These "tests" show if we can provide the arguments
    # to the workflow elements using Task and Param
    # classes.
    # Here we just read 2 types of data tables,
    # a number of variants each.
    test1 = Task$new(
        method = ResourceByEntity,
        name = NULL,
        res_by_entity_input_param
    ),
    
    test2 = Task$new(
        method = IntercellCategoriesPairwise,
        name = NULL,
        con_enrich_input_param
    ),
    
    # actual tasks of the workflow
    
    resource_coverage = Task$new(
        method = ResourceCoverage,
        name = NULL,
        con_enrich_input_param,
        only_main_classes = Param$new(TRUE, FALSE)
    ),
    
    intercell_subcls_intersect = Task$new(
        method = SubclassesIntersectionSeries,
        name = NULL,
        Param$new(NULL),
        Param$new(FALSE
            #TRUE
        )
    ),
    
    complex_overlaps = Task$new(
        method = UpsetComplexes,
        name = NULL,
        Param$new(NULL)
    ),
    
    connection_enrichment = Task$new(
        method = ConnectionEnrichment,
        name = NULL,
        input_param = con_enrich_input_param,
        directed = Param$new(TRUE, FALSE),
        sign = Param$new(NULL, 'stim', 'inh'),
        only_main_classes = Param$new(FALSE, TRUE),
        heatmap_variables = Param$new('enrich-count', 'enrich', 'count')
    ),
    
    intercell_graph = Task$new(
        method = ConnectionGraph,
        name = NULL,
        con_enrich_input_param
    ),
    
    intercell_class_sizes = Task$new(
        method = IntercellClassSizesSeries,
        name = NULL,
        input_param = con_enrich_input_param
    ),
    
    enzyme_substrate_modtype = Task$new(
        method = EnzymeSubstrateModtype,
        name = NULL
    ),
    
    enzyme_substrate_modtype_dot = Task$new(
        method = EnzymeSubstrateModtypeDot,
        name = NULL
    ),
    
    enzyme_substrate_self = Task$new(
        method = EnzymeSubstrateSelf,
        name = NULL,
        bar = Param$new(FALSE, TRUE),
        log_y = Param$new(FALSE, TRUE)
    ),
    
    enzyme_substrate_shared = Task$new(
        method = EnzymeSubstrateShared,
        name = NULL
    ),
    
    enzyme_substrate_numof_res = Task$new(
        method = EnzymeSubstrateNumofResources,
        name = NULL
    ),
    
    enzymes_per_substrate = Task$new(
        method = EnzymesPerSubstrate,
        name = NULL
    ),
    
    substrates_per_enzyme = Task$new(
        method = SubstratesPerEnzyme,
        name = NULL
    ),
    
    substrate_numof_sites = Task$new(
        method = SubstrateNumofSites,
        name = NULL
    ),
    
    complexes_by_resource = Task$new(
        method = ComplexesByResource,
        name = NULL,
        log_y = Param$new(TRUE, FALSE),
        bar = Param$new(TRUE, FALSE)
    ),
    
    complexes_numof_comp = Task$new(
        method = ComplexesNumofComponents,
        name = NULL
    ),
    
    refs_by_complex = Task$new(
        method = RefsPerComplex,
        name = NULL
    ),
    
    complexes_by_ref = Task$new(
        method = ComplexesPerRef,
        name = NULL
    ),
    
    complexes_by_components = Task$new(
        method = ComplexesByComponents,
        name = NULL
    ),
    
    network_sizes = Task$new(
        method = NetworkSizeDot,
        name = 'Network sizes',
        input_param = network_param
    ),
    
    network_directions = Task$new(
        method = NetworkDirectionsDot,
        name = 'Network directions',
        input_param = network_param
    ),
    
    network_coverage = Task$new(
        method = NetworkCoverageDot,
        name = 'Network coverages',
        input_param = network_param
    ),
    
    intercell_class_sizes_dot = Task$new(
        method = IntercellClassSizesDots,
        name = 'Intercell category sizes'
    )
    
)


omnipath2_run <- function(){
    
    for(task in omnipath2_workflow){
        
        task$run()
        
    }
    
}
