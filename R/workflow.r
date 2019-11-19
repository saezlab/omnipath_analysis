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

require(R6)


con_enrich_input_param <- ProductParam$new(
    list('curated', 'omnipath'),
    list('proteins'),
    list(
        'all-class-levels',
        'above_main-main-misc-small_main'
    )
)

res_by_entity_input_param <- Param$new('omnipath', 'curated')


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
        res_by_entity_input_param
    ),
    
    intercell_subcls_intersect = Task$new(
        method = SubclassesIntersectionSeries,
        name = NULL,
        Param$new(NULL),
        Param$new(TRUE, FALSE)
    ),
    
    complex_overlaps = Task$new(
        method = UpsetComplexes,
        name = NULL,
        Param$new(NULL)
    ),
    
    connection_enrichment = Task$new(
        method = ConnectionEnrichment,
        name = NULL,
        input_param = con_enrich_input_param
    ),
    
    intercell_graph = Task$new(
        method = ConnectionGraph,
        name = NULL,
        con_enrich_input_param
    ),
    
    category_sizes = Task$new(
        method = CategorySizesSeries,
        name = NULL,
        input_param = con_enrich_input_param
    )
    
)


omnipath2_run <- function(){
    
    for(task in omnipath2_workflow){
        
        task$run()
        
    }
    
}
