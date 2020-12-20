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
require(dplyr)
require(stringr)
require(R6)


NetworkCoverageDot <- R6::R6Class(
    
    'NetworkCoverageDot',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(input_param, only_totals = FALSE, ...){

            self$input_param <- input_param
            self$only_totals <- only_totals

            super$initialize(
                data = self$data,
                name = fig_network_coverage,
                fname_param = list(
                    input_param,
                    `if`(
                        self$only_totals,
                        'totals',
                        'by-resource'
                    )
                ),
                height = 9,
                width = 6
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            sc <- max(self$data$n_network / 100)
            
            self$plt <- ggplot(
                    self$data,
                    aes(
                        y = coverage_pct,
                        x = resource_label,
                        color = group,
                        shape = group
                    )
                ) +
                geom_col(
                    data = self$data %>%
                        group_by(resource_label, category) %>%
                        summarize_all(first),
                    mapping = aes(
                        y = n_network / sc,
                        x = resource_label,
                    ),
                    fill = '#CCCCCC',
                    color = NA
                ) +
                geom_point(
                    size = 2.3,
                    alpha = .8,
                    position = position_dodge(width = .4),
                    stroke = 0
                ) +
                scale_shape_manual(
                    name = 'Protein\ncategories',
                    values = omnipath2_settings$get('shapes')
                ) +
                scale_y_continuous(
                    sec.axis = sec_axis(
                        ~.*sc,
                        name = 'Number of proteins\nin resources'
                    )
                ) +
                scale_color_manual(
                    name = 'Protein\ncategories',
                    values = omnipath2_settings$get(palette2)
                ) +
                ylab('Coverage [%]') +
                xlab('Resources') +
                coord_flip()

                if(!self$only_totals){

                    self$plt <- self$plt +
                        facet_grid(
                            category~.,
                            scales = 'free_y',
                            space = 'free_y'
                        )

                }
            
            invisible(self)
            
        }
        
    ),
    
    private = list(
        
        load_data = function(){

            self$data <- `if`(
                self$only_totals,
                Networks$new(coverage = TRUE)$data,
                NetworkCoverage$new(input_param = self$input_param)$data
            )

            self$sizes <- `if`(
                self$only_totals,
                Networks$new()$data,
                NetworkSize$new(input_param = self$input_param)$data
            )

            invisible(self)

        },


        setup = function(){
            
            private$load_data()

            resource_order <- self$sizes %>%
                {`if`(
                    self$only_totals,
                        arrange(., dataset, desc(n)),
                        arrange(., category, desc(n))
                )} %>%
                pull(resource_label) %>%
                unique() %>%
                sub('miRNA-', 'miRNA ', .) %>%
                sub(' \\(total\\)', '', .)
            
            if(!self$only_totals){

                categories <- self$data %>%
                    filter(resource_type != 'resource') %>%
                    pull(resource_label) %>%
                    unique()

                resource_order <- c(
                    categories,
                    setdiff(resource_order, categories)
                )

            }

            self$data <- self$data %>%
                filter(resource_label != 'KEGG') %>%
                rename(category = data_model) %>%
                mutate(
                    resource_type = factor(
                        resource_type,
                        levels = c('total', 'data_model', 'resource'),
                        ordered = TRUE
                    ),
                    resource_label = factor(
                        resource_label,
                        levels = resource_order,
                        ordered = TRUE
                    ),
                    category = recode(category, All = 'OmniPath')
                ) %>%
                arrange(
                    category
                ) %>%
                mutate(coverage_pct = coverage / n_group * 100) %>%
                filter(!is.na(category) | resource_type != 'resource') %>%
                mutate(
                    category = ifelse(
                        is.na(category),
                        as.character(resource_label),
                        category
                    ),
                    category = factor(
                        category,
                        levels = unique(category),
                        ordered = TRUE
                    )
                ) %>%
                filter(!is.na(resource_label))

            self$height <- 1 + .2 * length(unique(self$data$resource_label))
            
            invisible(self)
            
        }
        
    )
    
)


NetworkSizeDot <- R6::R6Class(
    
    'NetworkSizeDot',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(input_param, only_totals = FALSE, ...){
            
            self$input_param <- input_param
            self$only_totals <- only_totals
            
            super$initialize(
                data = self$data,
                name = fig_network_size,
                fname_param = list(
                    input_param,
                    `if`(
                        self$only_totals,
                        'totals',
                        'by-resource'
                    )
                ),
                height = 9,
                width = 6
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            pal <- omnipath2_settings$get(palette2)
            
            StackedGroupedBarDot$new(
                obj = self,
                data = self$data,
                xvar = resource_label,
                yvar = n,
                fillvar = var,
                alphavar = shared,
                xlab = 'Resources',
                ylab = 'Number of components',
                log_y = TRUE,
                bar = FALSE,
                color_values = `names<-`(
                    pal,
                    c('entities', 'interactions_0')
                ),
                color_labels = c(
                    entities = 'Nodes',
                    interactions_0 = 'Interactions'
                ),
                legend_title = 'Network components',
                scale_alpha_param = list(
                    values = c(
                        `FALSE` = .8,
                        `TRUE` = .4
                    ),
                    labels = c(
                        `FALSE` = 'Total',
                        `TRUE` = 'Shared'
                    ),
                    guide = guide_legend(title = '')
                ),
                shape = 16
            )
            
            if(!self$only_totals){

                self$plt <- self$plt +
                    facet_grid(
                        category~.,
                        scales = 'free_y',
                        space = 'free_y'
                    )

            }
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        load_data = function(){

            self$data <- `if`(
                self$only_totals,
                Networks$new()$data,
                NetworkSize$new(input_param = self$input_param)$data
            )

            invisible(self)

        },


        setup = function(){

            private$load_data()

            categories <- self$data %>%
                filter(category == resource_label) %>%
                pull(resource_label) %>%
                unique()

            self$data <- self$data %>%
                {`if`(
                    self$only_totals,
                    arrange(., dataset, desc(n)),
                    arrange(., category, desc(n))
                )} %>%
                filter(resource_label != 'KEGG') %>%
                mutate(
                    resource_label = factor(
                        resource_label,
                        levels = c(
                            categories,
                            setdiff(unique(resource_label), categories)
                        ),
                        ordered = TRUE
                    )
                ) %>%
                filter(var %in% c('entities', 'interactions_0'))
            
            self$height <- 1 + .2 * length(unique(self$data$resource_label))
            
            invisible(self)
            
        }
        
    )
    
)


NetworkDirectionsDot <- R6::R6Class(
    
    'NetworkDirectionsDot',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(input_param, only_totals = FALSE, ...){
            
            self$input_param <- input_param
            self$only_totals <- only_totals
            
            self$vars <- c(
                'interactions_non_directed_0',
                'interactions_directed',
                'interactions_positive',
                'interactions_negative',
                'interactions_0'
            )
            
            super$initialize(
                data = self$data,
                name = fig_network_dir,
                fname_param = list(
                    input_param,
                    `if`(
                        self$only_totals,
                        'totals',
                        'by-resource'
                    )
                ),
                height = 9,
                width = 6
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            pal <- omnipath2_settings$get(palette2)
            
            StackedGroupedBarDot$new(
                obj = self,
                data = self$data,
                xvar = resource_label,
                yvar = n,
                fillvar = var,
                xlab = 'Resources',
                ylab = 'Number of interactions',
                log_y = TRUE,
                bar = FALSE,
                color_values = `names<-`(
                    pal,
                    self$vars
                ),
                color_labels = `names<-`(
                    c(
                        'Undirected',
                        'Directed',
                        'Positive',
                        'Negative',
                        'Total'
                    ),
                    self$vars
                ),
                legend_title = 'Interaction types',
                shape = -1,
                size = 2.5,
                position = position_dodge(width = .4)
            )
            
            if(!self$only_totals){

                self$plt <- self$plt +
                    facet_grid(
                        category~.,
                        scales = 'free_y',
                        space = 'free_y'
                    )

            }
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        load_data = function(){

            self$data <- `if`(
                self$only_totals,
                Networks$new()$data,
                NetworkSize$new(input_param = self$input_param)$data
            )

            invisible(self)

        },


        setup = function(){
            
            private$load_data()

            categories <- self$data %>%
                filter(category == resource_label) %>%
                pull(resource_label) %>%
                unique()

            self$data <- self$data %>%
            filter(resource_label != 'KEGG') %>%
            {`if`(
                self$only_totals,
                    arrange(., dataset, desc(n)),
                    arrange(., category, desc(n))
            )} %>%
            mutate(
                resource_label = factor(
                    resource_label,
                    levels = `if`(
                        self$only_totals,
                        unique(resource_label),
                        c(
                            categories,
                            setdiff(unique(resource_label), categories)
                        )
                    ),
                    ordered = TRUE
                ),
                var = factor(
                    var,
                    levels = rev(self$vars),
                    ordered = TRUE
                )
            ) %>%
            filter(var %in% self$vars & !shared) %>%
            filter(n > 0)
            
            self$height <- 1 + .2 * length(unique(self$data$resource))
            
            invisible(self)
            
        }
        
    )
    
)
