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
require(dplyr)
require(rlang)
require(R6)


ComplexBase <- R6::R6Class(
    
    'ComplexBase',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            name,
            expand_members = FALSE,
            xlab_vertical = TRUE,
            ...
        ){
            
            self$expand_members <- expand_members
            theme_args <- `if`(xlab_vertical, x_vertical_labels(), list())
            
            args <- modifyList(
                list(
                    data = self$data,
                    name = enquo(name),
                    width = 5,
                    theme_args = theme_args
                ),
                list(...)
            )
            
            do.call(super$initialize, args)
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            self$complexes <- Complexes$new(
                expand_members = self$expand_members,
                expand_resources = TRUE
            )$data
            
            invisible(self)
            
        }
        
    )
    
)


ComplexesByResource <- R6::R6Class(
    
    'ComplexesByResource',
    
    inherit = ComplexBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(name = fig_cplex_by_resource)
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(
                    self$complexes,
                    aes(y = n_complexes, x = resource)
                ) +
                geom_col(aes(fill = shared)) +
                scale_fill_manual(
                    values = c(
                        `TRUE` = '#4268B3',
                        `FALSE` = '#B3C5E9'
                    ),
                    labels = c(
                        `TRUE` = 'Shared',
                        `FALSE` = 'Unique'
                    ),
                    name = 'Number of\ncomplexes'
                ) +
                xlab('Resouces') +
                ylab('Complexes')
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$complexes <- self$complexes %>%
                group_by(complex_id) %>%
                mutate(shared = n() > 1) %>%
                ungroup()
            
            self$complexes <- bind_rows(
                
                self$complexes %>%
                    group_by(resource, shared) %>%
                    mutate(n_complexes = n()) %>%
                    summarize_all(first) %>%
                    ungroup(),
                
                self$complexes %>%
                    group_by(complex_id) %>%
                    summarize_all(first) %>%
                    ungroup() %>%
                    group_by(shared) %>%
                    mutate(
                        n_complexes = n(),
                        resource = 'Total'
                    ) %>%
                    summarize_all(first) %>%
                    ungroup()
                
            ) %>%
            group_by(resource) %>%
            mutate(n_complexes_total = sum(n_complexes)) %>%
            ungroup() %>%
            arrange(desc(n_complexes_total)) %>%
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