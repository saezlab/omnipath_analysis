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
                mutate(
                    resource = gsub('Signor', 'SIGNOR', resource)
                ) %>%
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


Networks <- R6::R6Class(

    'Networks',

    lock_objects = FALSE,

    public = list(

        initialize = function(
                networks = NULL,
                only_totals = TRUE,
                coverage = FALSE,
                ...
            ){

            self$networks <- `if`(
                is.null(networks),
                omnipath2_settings$get(network_datasets),
                networks
            )
            self$only_totals <- only_totals
            self$coverage <- coverage
            self$loader <- `if`(
                self$coverage,
                NetworkCoverage,
                NetworkSize
            )

            self$collect_networks()

            invisible(self)

        },


        collect_networks = function(){

            data <- NULL

            for(dataset in self$networks){

                reader <- self$loader$new(input_param = dataset)
                this_data <- reader$data %>%
                    mutate(dataset = dataset)
                data <- data %>% bind_rows(this_data)

            }

            self$data <- data %>%
                {`if`(
                    self$only_totals,
                    {`if`(
                        self$coverage,
                        filter(., resource_type != 'resource'),
                        filter(., grepl('total', resource))
                    )} %>%
                    mutate(
                        dataset = factor(
                            dataset,
                            levels = c(
                                'tf_mirna',
                                'mirna_mrna',
                                'tf_target',
                                'omnipath'
                            ),
                            ordered = TRUE
                        ),
                        resource_label = self$get_resource_label(resource),
                        dataset_label = map_chr(
                            as.character(dataset),
                            function(x){
                                omnipath2_settings$get(
                                    network_dataset_labels
                                )[[x]]
                            }
                        )
                    ) %>%
                    {`if`(
                        self$coverage,
                        mutate(
                            .,
                            resource_label = ifelse(
                                resource_type == 'total',
                                dataset_label,
                                resource_label
                            ),
                            resource_label = ifelse(
                                startsWith(resource_label, 'miRNA'),
                                resource_label,
                                str_to_sentence(resource_label)
                            )
                        ) %>%
                        filter(
                            resource != 'Activity flow' |
                            dataset == 'omnipath'
                        ),
                        filter(
                            .,
                            category != 'Activity flow' |
                            dataset == 'omnipath'
                        ) %>%
                        arrange(dataset, desc(category)) %>%
                        mutate(
                            resource = factor(
                                resource,
                                levels = unique(resource),
                                ordered = TRUE
                            )
                        )
                    )},
                    .
                )}

            invisible(self)

        },


        get_resource_label = function(label){

            label %>%
            str_replace(' total', '') %>%
            str_replace(
                'Post translational',
                'Post translational (total)'
            ) %>%
            str_replace('Mirna', 'miRNA') %>%
            str_replace(' trans', '-trans')

        }

    )

)


IntercellNetworkByResource <- R6::R6Class(

    'IntercellNetworkByResource',

    inherit = Reader,

    lock_objects = FALSE,

    public = list(

        initialize = function(...){

            super$initialize(input_intercell_network_by_resource_tsv)

            invisible(self)

        },


        preprocess = function(...){

            self$data <- self$data %>%
                select(id_a, id_b, database_a, sources) %>%
                group_by(database_a, id_a, id_b) %>%
                summarize_all(first) %>%
                group_by(id_a, id_b) %>%
                mutate(
                    con_is_unique = n_distinct(database_a) <= 2
                ) %>%
                ungroup() %>%
                group_by(id_a) %>%
                mutate(
                    a_is_unique = n_distinct(database_a) <= 2
                ) %>%
                ungroup() %>%
                group_by(id_b) %>%
                mutate(
                    b_is_unique = n_distinct(database_a) <= 2
                ) %>%
                ungroup() %>%
                group_by(database_a, a_is_unique) %>%
                mutate(
                    n_a_unique = n_distinct(id_a) * a_is_unique
                ) %>%
                ungroup() %>%
                group_by(database_a, b_is_unique) %>%
                mutate(
                    n_b_unique = n_distinct(id_b) * b_is_unique
                ) %>%
                ungroup() %>%
                group_by(database_a) %>%
                mutate(
                    n_con_unique = sum(con_is_unique),
                    n_a_unique = max(n_a_unique),
                    n_a = n_distinct(id_a),
                    n_b_unique = max(n_a_unique),
                    n_b = n_distinct(id_b),
                    n_con = n_distinct(id_a, id_b)
                ) %>%
                summarize_all(first) %>%
                select(
                    resource = database_a,
                    n_con,
                    n_con_unique,
                    n_a,
                    n_a_unique,
                    n_b,
                    n_b_unique
                ) %>%
                arrange(desc(n_con)) %>%
                mutate(
                    resource = factor(
                        resource,
                        levels = unique(resource),
                        ordered = TRUE
                    )
                ) %>%
                gather(key = 'var', value = 'cnt', -resource) %>%
                mutate(
                    uni = grepl('unique', var, fixed = TRUE),
                    obj = ifelse(
                        grepl('con', var, fixed = TRUE),
                        'con',
                        ifelse(
                            grepl('n_a', var, fixed = TRUE),
                            'transmitter',
                            'receiver'
                        )
                    )
                ) %>%
                mutate(
                    obj = factor(
                        obj,
                        levels = c('con', 'transmitter', 'receiver'),
                        ordered = TRUE
                    )
                )

            invisible(self)

        }

    )

)
