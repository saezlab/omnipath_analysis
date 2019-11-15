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


list_product <- function(lst){
    
    unlist(
        apply(
            unname(expand.grid(lst)),
            1,
            list
        ),
        recursive = FALSE
    )
    
}


ProductParam <- R6::R6Class(
    
    'ProductParam',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            self$param <- unlist(
                apply(
                    expand.grid(list(...)),
                    1,
                    list
                ),
                recursive = FALSE
            )
            
            invisible(self)
            
        }
        
    )
    
)


Param <- R6::R6Class(
    
    'Param',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            self$param <- list(list(...))
            
        }
        
    )
    
)


Task <- R6::R6Class(
    
    'Task',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            method,
            name = NULL,
            ...
        ){
        
            self$method <- method
            self$param <- list(...)
            self$name <- `if`(
                is.null(name),
                as.character(expr(enquo(method))),
                name
            )
            
            private$set_param()
            
            invisible(self)
            
        },
        
        
        run = function(){
            
            for(param in private$product_param()){
                
                self$current_param <- param
                private$run_one()
                
            }
            
            invisible(self)
            
        }
        
    ),
    
    private = list(
        
        set_param = function(){
            
            self$param <- case_when(
                length(self$param) == 0 ~ list(Param$new()),
                is.list(self$param) ~ self$param,
                TRUE ~ list(self$param)
            )
            
            invisible(self)
            
        },
        
        
        product_param = function(){
            
            do.call(
                list_product,
                lapply(
                    self$param,
                    function(p){p$param}
                )
            )
        },
        
        
        run_one = function(){
            
            self$worker <- do.call(
                self$method$new,
                self$current_param
            )
            
            invisible(self)
            
        }
        
    )
    
)
