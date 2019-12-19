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
require(grid)
require(gridExtra)
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
            width = 5,
            ...
        ){
            
            self$expand_members <- expand_members
            self$expand_resources <- expand_resources
            theme_args <- `if`(xlab_vertical, x_vertical_labels(), list())
            args <- list(...)
            theme_args <- modifyList(
                theme_args,
                `if`(
                    !is.null(args$theme_args),
                    args$theme_args,
                    list()
                )
            )
            
            args <- modifyList(
                list(
                    data = self$data,
                    name = enquo(name),
                    width = width,
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
        
        initialize = function(log_y = FALSE, bar = TRUE){
            
            self$log_y <- log_y
            self$bar <- bar
            
            super$initialize(
                name = fig_cplex_by_resource,
                fname_param = list(
                    `if`(log_y, 'log', 'stacked'),
                    `if`(bar, 'bar', 'dot')
                ),
                height = 4,
                width = 4,
                xlab_vertical = FALSE
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            StackedGroupedBarDot$new(
                obj = self,
                data = self$complexes,
                xvar = resource,
                yvar = n_complexes,
                fillvar = shared,
                xlab = 'Resource',
                ylab = sprintf(
                    'Complexes%s',
                    `if`(self$log_y, ' (log)', '')
                ),
                log_y = self$log_y,
                bar = self$bar,
                color_values = `names<-`(
                    omnipath2_settings$get(two_shades_1),
                    c(TRUE, FALSE)
                ),
                color_labels = c(
                    `TRUE` = 'Shared',
                    `FALSE` = 'Unique'
                ),
                legend_title = 'Number of\ncomplexes'
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            secondary_sources = c('NetworkBlast', 'CFinder')
            
            super$setup()
            
            self$complexes <-self$complexes %>%
                group_by(complex_id) %>%
                mutate(
                    shared = length(setdiff(resource, secondary_sources)) > 1
                ) %>%
                ungroup() %>%
                mutate(resource = as.character(resource)) %>%
                full_join(
                    expand.grid(
                        resource =unique(self$complexes$resource),
                        shared = c(TRUE, FALSE)
                    ) %>%
                    mutate(resource = as.character(resource)),
                    by = c('resource', 'shared')
                )
            
            self$complexes <- bind_rows(
                
                self$complexes %>%
                    group_by(resource, shared) %>%
                    mutate(n_complexes = n_distinct(complex_id)) %>%
                    summarize_all(first) %>%
                    ungroup(),
                
                self$complexes %>%
                    group_by(complex_id) %>%
                    summarize_all(first) %>%
                    ungroup() %>%
                    group_by(shared) %>%
                    mutate(
                        n_complexes = n(),
                        resource = 'OmniPath'
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
                xlab_vertical = FALSE,
                width = 7
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            DistDensHist$new(
                obj = self,
                data = self$complexes,
                x = n_members,
                ylab = 'Complexes',
                xlab = 'Number of components\nin the complex',
                density_adjust = 3
            )
            
            invisible(self)
            
        },
            
        save = function(){
            
            super$save(print_open = FALSE)
            
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
                xlab_vertical = FALSE,
                width = 7
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            DistDensHist$new(
                obj = self,
                data = self$complexes,
                x = n_complexes,
                ylab = 'Components',
                xlab = 'Number of complexes\nincluding the component',
                density_adjust = 3.5
            )
            
            invisible(self)
            
        },
        
        save = function(){
            
            super$save(print_open = FALSE)
            
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


ComplexesPerRef <- R6::R6Class(
    
    'ComplexesPerRef',
    
    inherit = ComplexBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_complexes_by_ref,
                expand_resources = FALSE,
                xlab_vertical = FALSE,
                width = 7
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            DistDensHist$new(
                obj = self,
                data = self$complexes,
                x = n_complexes,
                ylab = 'References',
                xlab = 'Number of complexes\nfrom a reference',
                density_adjust = 4
            )
            
            invisible(self)
            
        },
        
        save = function(){
            
            super$save(print_open = FALSE)
            
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


RefsPerComplex <- R6::R6Class(
    
    'RefsPerComplex',
    
    inherit = ComplexBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_refs_by_cplex,
                expand_resources = FALSE,
                xlab_vertical = FALSE,
                width = 7
            )
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            DistDensHist$new(
                obj = self,
                data = self$complexes,
                x = n_refs,
                ylab = 'Complexes',
                xlab = 'Number of references\nfor a complex',
                density_adjust = 5
            )
            
            invisible(self)
            
        },
        
        save = function(){
            
            super$save(print_open = FALSE)
            
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
