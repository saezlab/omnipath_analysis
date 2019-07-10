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
require(R6)


PlotSeries <- R6::R6Class(
    
    'PlotSeries',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            data,
            slice_var,
            plotter,
            name,
            plot_args = list(),
            label_mapping = NULL,
            exclude_levels = NULL,
            label_levels = NULL
        ){
            
            self$data <- data
            self$name <- name
            self$slice_var <- enquo(slice_var)
            self$exclude_levels <- exclude_levels
            self$label_levels <- label_levels
            self$plotter <- plotter
            self$plot_args <- plot_args
            
            self$main()
            
            invisible(self)
            
        },
        
        
        main = function(){
            
            self$preprocess()
            self$set_levels()
            
            #self$plot()
            
            invisible(self)
            
        },
        
        
        preprocess = function(...){
            
            invisible(self)
            
        },
        
        
        iter_levels = function(){
            
            idx <- 0
            
            function(reset = FALSE){
                
                idx <<- `if`(reset, 1, idx + 1)
                
                `if`(
                    idx <= length(self$use_levels),
                    self$use_levels[idx],
                    NULL
                )
                
            }
            
        },
        
        
        iter_slices = function(){
            
            level_iterator <- self$iter_levels()
            
            function(reset = FALSE){
                
                level <- level_iterator(reset = reset)
                
                if(!is.null(level)){
                    
                    self$level <- level
                    self$slice_label <- `if`(
                        level %in% self$label_levels,
                        self$label_levels[[level]],
                        level
                    )
                    self$slice_name <- sprintf(
                        '%s-by-%s__%s',
                        self$name,
                        quo_name(self$slice_var),
                        level
                    )
                    self$slice <- self$data %>%
                        filter(!!self$slice_var == level)
                    
                }else{
                    
                    self$level <- NULL
                    self$slice_label <- NULL
                    self$slice_name <- NULL
                    self$slice <- NULL
                    
                }
                
            }
            
        },
        
        
        plot = function(){
            
            slice_iterator <- self$iter_slices()
            
            slice_iterator()
            
            while(!is.null(self$slice)){
                
                private$single_plot()
                
                slice_iterator()
                
            }
            
            invisible(self)
            
        },
        
        
        set_levels = function(){
            
            self$use_levels <- setdiff(
                (
                    self$data %>%
                    select(!!self$slice_var) %>%
                    c %>% unique %>% unlist
                ),
                self$exclude_levels
            )
            
        }
        
    ),
    
    private = list(
        
        single_plot = function(){
            
            do.call(
                self$plotter$new,
                c(
                    list(
                        data = self$slice,
                        name = self$slice_name
                    ),
                    self$plot_args
                )
            )
            
        }
        
    )
    
)


CategoriesPairwisePlotSeries <- R6::R6Class(
    
    'CategoriesPairwisePlotSeries',
    
    inherit = PlotSeries,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(slice_var, plotter, name, data = NULL, ...){
            
            slice_var <- enquo(slice_var)
            
            private$ensure_data(data)
            
            super$initialize(
                data = self$data,
                slice_var = !!slice_var,
                plotter = plotter,
                name = name,
                ...
            )
            
        }
        
    ),
    
    
    private = list(
        
        ensure_data = function(data){
            
            self$data <- `if`(
                is.null(data),
                IntercellCategoriesPairwise$new()$data,
                data
            )
            
            invisible(self)
            
        }
        
    )
    
)
