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
            expand_resources = TRUE,
            xlab_vertical = TRUE,
            ...
        ){
            
            self$expand_members <- expand_members
            self$expand_resources <- expand_resources
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
                expand_resources = self$expand_resources,
                keep_references = TRUE,
                keep_stoichiometry = TRUE
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
            
            secondary_sources = c('NetworkBlast', 'CFinder')
            
            super$setup()
            
            self$complexes <- self$complexes %>%
                group_by(complex_id) %>%
                mutate(
                    shared = length(setdiff(resource, secondary_sources)) > 1
                ) %>%
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


ComplexesNumofComponents <- R6::R6Class(
    
    'ComplexesNumofComponents',
    
    inherit = ComplexBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_cplex_n_comp,
                expand_resources = FALSE,
                xlab_vertical = FALSE
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(self$complexes) +
                stat_bin(
                        aes(x = n_members, y = cumsum(..count..)),
                        geom = 'step',
                        binwidth = .01
                ) +
                xlab('Number of components') +
                ylab('Complexes') +
                scale_x_log10() +
                annotation_logticks(sides = 'b')
            
            invisible(self)
            
        }
        
    )
    
)


ComplexesByComponents <- R6::R6Class(
    
    'ComplexesByComponents',
    
    inherit = ComplexBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_cplex_by_comp,
                expand_resources = FALSE,
                expand_members = TRUE,
                xlab_vertical = FALSE
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(self$complexes) +
                stat_bin(
                        aes(x = n_members, y = cumsum(..count..)),
                        geom = 'step',
                        binwidth = .01
                ) +
                xlab('Number of complexes') +
                ylab('Components') +
                scale_x_log10() +
                annotation_logticks(sides = 'b')
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$complexes <- self$complexes %>%
                group_by(members) %>%
                mutate(n_complexes = n_distinct(complex_id)) %>%
                summarize_all(first) %>%
                ungroup()
            
            invisible(self)
            
        }
        
    )
    
)


RefsPerComplex <- R6::R6Class(
    
    'RefsPerComplex',
    
    inherit = ComplexBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_refs_by_cplex,
                expand_resources = FALSE,
                xlab_vertical = FALSE
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(self$complexes) +
                stat_bin(
                        aes(x = n_complexes, y = cumsum(..count..)),
                        geom = 'step',
                        binwidth = .01
                ) +
                xlab('Number of references') +
                ylab('Complexes') +
                scale_x_log10() +
                annotation_logticks(sides = 'b')
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$complexes <- self$complexes %>%
                separate_rows(references, sep = ';') %>%
                group_by(references) %>%
                mutate(n_complexes = n_distinct(complex_id)) %>%
                summarize_all(first) %>%
                ungroup()
            
            invisible(self)
            
        }
        
    )
    
)


ComplexesPerRef <- R6::R6Class(
    
    'ComplexesPerRef',
    
    inherit = ComplexBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_complexes_by_ref,
                expand_resources = FALSE,
                xlab_vertical = FALSE
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(self$complexes) +
                stat_bin(
                        aes(x = n_refs, y = cumsum(..count..)),
                        geom = 'step',
                        binwidth = .01
                ) +
                xlab('Number of complexes') +
                ylab('References') +
                scale_x_log10() +
                annotation_logticks(sides = 'b')
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$complexes <- self$complexes %>%
                separate_rows(references, sep = ';') %>%
                group_by(complex_id) %>%
                mutate(n_refs = n_distinct(references)) %>%
                summarize_all(first) %>%
                ungroup()
            
            invisible(self)
            
        }
        
    )
    
)
