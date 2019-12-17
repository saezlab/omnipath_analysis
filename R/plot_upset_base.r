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
require(UpSetR)
require(grid)
require(R6)


UpsetBase <- R6::R6Class(
    
    'UpsetBase',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
                name,
                cat_col,
                ent_col,
                fname_param = list(),
                #title = NULL,
                data = NULL,
                ...
            ){
            
            self$data <- `if`(is.null(self$data), data, self$data)
            self$cat_col <- enquo(cat_col)
            self$ent_col <- enquo(ent_col)
            #self$title <- title
            
            super$initialize(
                data = data,
                name = UQ(name),
                fname_param = fname_param,
                ...
            )
            
            invisible(self)
            
        },
        
        
        plot = function(...){
            
            self$plt <- do.call(
                upset,
                c(
                    list(
                        fromList(self$categories)
                    ),
                    list(...)
                ),
                
            )
            
            invisible(self)
            
        },
        
        
        post_plot = function(){
            
            if(!is.null(self$title)){
                
                grid.text(
                    self$title,
                    x = 0.65,
                    y = 0.95,
                    gp = gpar(fontsize = 20)
                )
                
            }
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            private$ensure_data()
            private$collect_categories()
            
        },
        
        
        ensure_data = function(){
            
            if(is.null(self$data)){
                
                stop('omnipath2::UpsetBase: No data provided.')
                
            }
            
        },
        
        
        collect_categories = function(){
            
            self$names <- self$data %>%
                select(!!self$cat_col) %>%
                unlist() %>%
                unique()
            
            self$categories <- setNames(
                sapply(
                    self$names,
                    function(name){
                        self$data %>%
                            filter(!!self$cat_col == name) %>%
                            select(!!self$ent_col) %>%
                            unlist() %>%
                            unique()
                    }
                ),
                self$names
            )
            
            if(is.null(self$plot_args$nsets)){
                
                self$plot_args$nsets <- length(self$names)
                
            }
            
        }
        
    )
    
)
