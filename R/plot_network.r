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
        
        initialize = function(input_param, ...){
            
            self$input_param <- input_param
            
            super$initialize(
                data = self$data,
                name = fig_network_coverage,
                fname_param = list(input_param),
                height = 9,
                width = 6
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            sc <- max(self$data$n_network / 100)
            
            self$plt <- ggplot(
                    self$data,
                    aes(y = coverage_pct, x = resource, color = group)
                ) +
                geom_col(
                    data = self$data %>%
                    group_by(resource) %>%
                    summarize_all(first),
                    mapping = aes(
                        y = n_network / sc,
                        x = resource,
                    ),
                    fill = '#CCCCCC',
                    color = NA
                ) +
                geom_point(
                    size = 2.3,
                    alpha = .8,
                    position = position_dodge(width = .4)
                ) +
                scale_y_continuous(
                    sec.axis = sec_axis(
                        ~.*sc,
                        name = 'Number of proteins\nin resources'
                    )
                ) +
                scale_color_manual(
                    guide = guide_legend(
                        title = 'Protein\ncategories'
                    ),
                    values = omnipath2_settings$get(palette2)
                ) +
                ylab('Coverage [%]') +
                xlab('Resources') +
                coord_flip()
            
            invisible(self)
            
        }
        
    ),
    
    private = list(
        
        setup = function(){
            
            self$data <- NetworkCoverage$new(
                input_param = self$input_param
            )$data %>%
            mutate(
                resource_type = factor(
                    resource_type,
                    levels = c('total', 'data_model', 'resource'),
                    ordered = TRUE
                )
            ) %>%
            arrange(resource_type, desc(n_network)) %>%
            mutate(coverage_pct = coverage / n_group * 100) %>%
            mutate(
                resource = factor(
                    resource,
                    levels = unique(resource),
                    ordered = TRUE
                )
            )
            
            invisible(self)
            
        }
        
    )
    
)


NetworkSizeDot <- R6::R6Class(
    
    'NetworkSizeDot',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(input_param, ...){
            
            self$input_param <- input_param
            
            super$initialize(
                data = self$data,
                name = fig_network_size,
                fname_param = list(input_param),
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
                xvar = resource,
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
                shape = 19
            )
            
            self$plt <- self$plt +
                facet_grid(
                    category~.,
                    scales = 'free_y',
                    space = 'free_y'
                )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            self$data <- NetworkSize$new(
                input_param = self$input_param
            )$data %>%
            filter(var %in% c('entities', 'interactions_0')) %>%
            arrange(desc(n)) %>%
            mutate(
                resource = factor(
                    resource,
                    levels = unique(resource),
                    ordered = TRUE
                )
            )
            
            self$height <- 1 + .2 * length(unique(self$data$resource))
            
            invisible(self)
            
        }
        
    )
    
)


NetworkDirectionsDot <- R6::R6Class(
    
    'NetworkDirectionsDot',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
                input_param, ...){
            
            self$input_param <- input_param
            
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
                fname_param = list(input_param),
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
                xvar = resource,
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
                shape = 16,
                size = 2.5,
                position = position_dodge(width = .4)
            )
            
            self$plt <- self$plt +
                facet_grid(
                    category~.,
                    scales = 'free_y',
                    space = 'free_y'
                )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            self$data <- NetworkSize$new(
                input_param = self$input_param
            )$data %>%
            filter(var %in% self$vars & !shared) %>%
            arrange(desc(n)) %>%
            mutate(
                resource = factor(
                    resource,
                    levels = unique(resource),
                    ordered = TRUE
                ),
                var = factor(
                    var,
                    levels = rev(self$vars),
                    ordered = TRUE
                )
            ) %>%
            filter(n > 0)
            
            self$height <- 1 + .2 * length(unique(self$data$resource))
            
            invisible(self)
            
        }
        
    )
    
)
