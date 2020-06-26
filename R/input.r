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
        
        initialize = function(
            input_param,
            entity_type = NULL,
            only_main_classes = FALSE,
            ...
        ){
            
            self$entity_type <- entity_type
            self$only_main_classes <- only_main_classes
            
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
            
            main_classes <- omnipath2_settings$get(intercell_main_classes)

            self$data <- self$data %>%
                {`if`(
                    is.null(self$entity_type) || !('entity' %in% names(.)),
                    .,
                    filter(., entity == self$entity_type)
                )} %>%
                mutate(
                    network_covers_cls0 = size0 / in_network_cls0 * 100,
                    network_covers_cls1 = size1 / in_network_cls1 * 100,
                    pct_of_parent_cls0 = size0 / size_parent0 * 100,
                    pct_of_parent_cls1 = size1 / size_parent1 * 100
                ) %>%
                arrange(desc(size0)) %>%
                mutate(
                    name0 = factor(
                        name0,
                        levels = unique(name0),
                        ordered = TRUE
                    )
                ) %>%
                {`if`(
                    self$only_main_classes,
                    filter(
                        .,
                        name0 %in% main_classes &
                        name1 %in% main_classes
                    ),
                    .
                )}
            
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


Complexes <- R6::R6Class(
    
    'Complexes',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            expand_members = FALSE,
            keep_references = TRUE,
            expand_resources = FALSE,
            keep_stoichiometry = FALSE,
            ...
        ){
            
            self$expand_members <- expand_members
            self$keep_references <- keep_references
            self$expand_resources <- expand_resources
            self$keep_stoichiometry <- keep_stoichiometry
            
            super$initialize(input_complexes_by_resource_tsv)
            
            invisible(self)
            
        },
        
        
        preprocess = function(){
            
            self$data <- self$data %>%
                {`if`(
                    self$keep_references,
                    .,
                    select(., -references)
                )} %>%
                {`if`(
                    self$keep_stoichiometry,
                    .,
                    select(., -stoichiometry)
                )} %>%
                mutate(members = sub('COMPLEX:', '', complex_id)) %>%
                separate_rows(members, sep = '_') %>%
                group_by(complex_id) %>%
                mutate(n_members = n_distinct(members)) %>%
                ungroup() %>%
                {`if`(
                    self$expand_members,
                    .,
                    group_by(., complex_id) %>%
                    mutate(
                        resource = paste(
                            unique(resource),
                            collapse = ';'
                        )
                    ) %>%
                    summarize_all(first) %>%
                    ungroup() %>%
                    select(-members)
                )} %>%
                {`if`(
                    self$expand_resources,
                    separate_rows(., resource, sep = ';'),
                    .
                )}
            
            invisible(self)
            
        }
        
    )
    
)


Annotations <- R6::R6Class(
    
    'Annotations',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            super$initialize(input_annotation_by_entity_tsv)
            
            invisible(self)
            
        }
        
    )
    
)


NetworkCoverage <- R6::R6Class(
    
    'NetworkCoverage',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(input_param, ...){
            
            do.call(
                super$initialize,
                c(
                    list(
                        name = quote(input_network_coverages_tsv)
                    ),
                    input_param
                )
            )
            
            invisible(self)
            
        }
        
    )
    
)


NetworkSize <- R6::R6Class(
    
    'NetworkSize',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(input_param, ...){
            
            vars <- c(
                'entities',
                'interactions_0',
                'interactions_non_directed_0',
                'interactions_directed',
                'interactions_positive',
                'interactions_negative'
            )
            
            do.call(
                super$initialize,
                c(
                    list(
                        name = quote(input_network_size_tsv)
                    ),
                    input_param
                )
            )
            
            resources <- rev(self$data %>% pull(resource))
            
            categories <- NULL
            category <- NA
            
            for(res in resources){
                
                if(endsWith(res, 'total')){
                    category <- sub(' total', '', res)
                }
                
                categories <- c(categories, category)
                
            }
            
            self$data <- self$data %>%
                add_column(category = rev(categories))
            
            data <- NULL
            
            for(var in vars){
                
                this_data <- self$data %>%
                    select(
                        resource = resource,
                        category = category,
                        n_total = sprintf(
                                '%s_n',
                                var
                        ),
                        n_shared = sprintf(
                                '%s_shared_within_database_category_n',
                                var
                        )
                    ) %>%
                    mutate(var = var)
                
                this_data <- bind_rows(
                    this_data %>%
                        select(
                            resource,
                            category,
                            n = n_total,
                            var
                        ) %>%
                        mutate(shared = FALSE),
                    this_data %>%
                        select(
                            resource,
                            category,
                            n = n_shared,
                            var
                        ) %>%
                        mutate(shared = TRUE)
                )
                
                data <- bind_rows(data, this_data)
                
            }
            
            self$data <- data
            
            invisible(self)
            
        }
        
    )
    
)