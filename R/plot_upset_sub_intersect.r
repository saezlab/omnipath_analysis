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

require(dplyr)
require(R6)


SubclassesIntersection <- R6::R6Class(
    
    'SubclassesIntersection',
    
    inherit = UpsetBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            data = NULL,
            name = fig_subcls_intersect,
            fname_param = list(),
            width = NULL,
            height = NULL,
            plot_args = list(),
            ...
        ){
            
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
                data = data,
                fname_param = fname_param,
                ent_col = entity_id,
                cat_col = resource_label,
                name = enquo(name),
                plot_args = plot_args,
                width = width,
                height = height,
                ...
            )
            
            invisible(self)
            
        }
        
    )
    
)


SubclassesIntersectionSeries <- R6::R6Class(
    
    'SubclassesIntersectionSeries',
    
    inherit = PlotSeries,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            data = NULL,
            complexes = FALSE,
            ...
        ){
            
            self$complexes <- complexes
            private$ensure_data(data)
            
            super$initialize(
                data = self$data,
                slice_var = parent,
                plotter = SubclassesIntersection,
                name = fig_subcls_intersect,
                fname_param = list(
                    `if`(self$complexes, 'complex', 'protein')
                ),
                width_by = resource_label,
                width_min = 4,
                width_step = .5
            )
            
            invisible(self)
            
        },
        
        preprocess = function(){
            
            self$data <- self$data %>%
                filter(
                    resource_label != '' &
                    class_type == 'sub'
                ) %>%
                {`if`(
                    self$complexes,
                    filter(., is_complex),
                    filter(., !is_complex)
                )}
            
        }
        
    ),
    
    
    private = list(
        
        ensure_data = function(data){
            
            self$data <- `if`(
                is.null(data),
                IntercellAnnotationByEntity$new()$data,
                data
            )
            
        }
        
    )
    
)
