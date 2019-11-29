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

require(rlang)
require(readr)
require(dplyr)
require(R6)


Reader <- R6::R6Class(
    
    'Reader',
    
    inherit = InputPath,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            name,
            ...
        ){
            
            super$initialize(
                name = UQ(enquo(name)),
                ...
            )
            
            self$main()
            
            invisible(self)
            
        },
        
        
        main = function(){
            
            self$read()
            self$preprocess()
            
            invisible(self)
            
        },
        
        
        read = function(){
            
            self$data <- suppressMessages(
                readr::read_tsv(
                    self$path,
                    progress = FALSE
                )
            )
            
            invisible(self)
            
        },
        
        
        preprocess = function(...){
            
            invisible(self)
            
        }
        
    )
    
)


IntercellCategoriesPairwise <- R6::R6Class(
    
    'IntercellCategoriesPairwise',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(input_param, entity_type = NULL, ...){
            
            self$entity_type <- entity_type
            
            do.call(
                super$initialize,
                c(
                    list(
                        name = quote(input_intercell_cat_pairwise_tsv)
                    ),
                    input_param
                )
            )
            
            invisible(self)
            
        },
        
        
        #' Here we add some more columns which are easy to calculate
        #' over the entire data frame. Also we sort the data frame by
        #' the size of the first category.
        preprocess = function(...){
            
            self$data <- self$data %>%
                {`if`(
                    is.null(self$entity_type) || !('entity' %in% names(.)),
                    .,
                    filter(., entity == self$entity_type)
                )} %>%
                mutate(
                    network_covers_cls0 = size_cls0 / in_network_cls0 * 100,
                    network_covers_cls1 = size_cls1 / in_network_cls1 * 100,
                    pct_of_parent_cls0 = size_cls0 / size_parent0 * 100,
                    pct_of_parent_cls1 = size_cls1 / size_parent1 * 100
                ) %>%
                arrange(desc(size_cls0)) %>%
                mutate(
                    name_cls0 = factor(
                        name_cls0,
                        levels = unique(name_cls0),
                        ordered = TRUE
                    )
                )
            
        }
        
    )
    
)


IntercellAnnotationByEntity <- R6::R6Class(
    
    'IntercellAnnotationByEntity',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            super$initialize(input_intercell_annotation_by_entity_tsv)
            
            invisible(self)
            
        },
        
        
        preprocess = function(...){
            
            self$data <- self$data %>%
                mutate(
                    is_complex = is_complex == 'True'
                )
            
            invisible(self)
            
        }
        
    )
    
)


AnnotationByEntity <- R6::R6Class(
    
    'IntercellAnnotationByEntity',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            super$initialize(input_intercell_annotation_by_entity_tsv)
            
            invisible(self)
            
        }
        
    )
    
)


EnzymeSubstrate <- R6::R6Class(
    
    'EnzymeSubstrate',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            super$initialize(input_enzyme_substrate_tsv)
            
            invisible(self)
            
        }
        
    )
    
)


ResourceByEntity <- R6::R6Class(
    
    'ResourceByEntity',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(input_param, ...){
            
            do.call(
                super$initialize,
                c(
                    list(
                        name = expr(input_resource_by_entity_tsv)
                    ),
                    input_param
                )
            )
            
            invisible(self)
            
        },
        
        
        preprocess = function(...){
            
            self$data <- self$data %>%
                mutate(
                    is_complex = is_complex == 'True'
                )
            
            invisible(self)
            
        }
        
    )
    
)


ComplexesByResource <- R6::R6Class(
    
    'ComplexesByResource',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(expand_members = FALSE, ...){
            
            self$expand_members <- expand_members
            
            super$initialize(input_complexes_by_resource_tsv)
            
            invisible(self)
            
        },
        
        
        preprocess = function(){
            
            self$data <- self$data %>%
                mutate(members = sub('COMPLEX:', '', complex_id)) %>%
                separate_rows(members, sep = '_') %>%
                group_by(complex_id) %>%
                mutate(n_members = n()) %>%
                ungroup() %>%
                {`if`(
                    self$expand_members,
                    .,
                    group_by(., complex_id) %>%
                    summarize_all(first) %>%
                    ungroup() %>%
                    select(-members)
                )}
            
            invisible(self)
            
        }
        
    )
    
)
