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


UpsetComplexes <- R6::R6Class(
    
    'UpsetComplexes',
    
    inherit = UpsetBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
                data = NULL,
                plot_args = list(),
                ...
            ){
            
            self$data <- data
            private$ensure_data()
            
            plot_args <- modifyList(
                list(
                    #title = 'Complexes by resource',
                    order.by = 'freq',
                    show.numbers = 'no',
                    mb.ratio = c(.55, .45),
                    scale.intersections = 'log10',
                    text.scale = 1.5
                ),
                plot_args
            )
            
            super$initialize(
                ent_col = complex_id,
                cat_col = resource,
                name = quote(fig_cplex_overlap),
                plot_args = plot_args,
                width = 7,
                height = 4,
                ...
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        ensure_data = function(){
            
            if(is.null(self$data)){
                
                private$read_data()
                
            }
            
        },
        
        
        read_data = function(){
            
            self$data <- Complexes$new()$data
            
        }
        
    )
    
)
