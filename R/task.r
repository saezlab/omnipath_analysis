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

require(R6)


list_product <- function(lst){
    
    lst_names <- names(lst)
    
    prod_lst <- as.list(
        unlist(
            apply(
                unname(expand.grid(lst)),
                1,
                list
            ),
            recursive = FALSE
        )
    )
    
    if(!is.null(lst_names)){
        
        prod_lst <- lapply(
            prod_lst,
            function(li){`names<-`(li, lst_names)}
        )
        
    }
    
    return(prod_lst)
    
}


ProductParam <- R6::R6Class(
    
    'ProductParam',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            self$param <- list_product(list(...))
            
            invisible(self)
            
        }
        
    )
    
)


Param <- R6::R6Class(
    
    'Param',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            self$param <- list(...)
            
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
                `if`(
                    class(method) == 'R6ClassGenerator',
                    method$classname,
                    as.character(expr(method))
                ),
                name
            )
            
            private$set_param()
            
            invisible(self)
            
        },
        
        
        run = function(){
            
            op2log(
                sprintf('Running task `%s`', self$name),
                label = 'workflow',
                0
            )
            
            for(param in private$product_param()){
                
                self$current_param <- param
                
                private$run_one()
                
            }
            
            invisible(self)
            
        }
        
    ),
    
    private = list(
        
        set_param = function(){
            
            self$param <- `if`(
                length(self$param) == 0,
                list(Param$new()),
                `if`(
                    is.list(self$param),
                    self$param,
                    list(self$param)
                )
            )
            
            invisible(self)
            
        },
        
        
        product_param = function(){
            
            list_product(
                sapply(
                    self$param,
                    function(p){p$param},
                    simplify = FALSE,
                    USE.NAMES = TRUE
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
