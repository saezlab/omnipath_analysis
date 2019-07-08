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

require(readr)
require(dplyr)
require(R6)


Reader <- R6::R6Class(
    
    'Reader',
    
    inherit = InputPath,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            name
        ){
            
            super$initialize(name = name)
            
            self$main()
            
            invisible(self)
            
        },
        
        
        main = function(){
            
            self$read()
            self$preprocess()
            
            invisible(self)
            
        },
        
        
        read = function(){
            
            self$data <- suppressMessages(readr::read_tsv(self$path))
            
            invisible(self)
            
        },
        
        
        #' Here we add some more columns which are easy to calculate
        #' over the entire data frame. Also we sort the data frame by
        #' the size of the first category.
        preprocess = function(...){
            
            self$data <- self$data %>%
                mutate(
                    omnipath_covers_cls0 = size_cls0 / omnipath0 * 100,
                    omnipath_covers_cls1 = size_cls1 / omnipath1 * 100,
                    pct_of_parent_cls0 = size_cls0 / size_parent0 * 100,
                    pct_of_parent_cls1 = size_cls1 / size_parent1 * 100
                ) %>%
                arrange(desc(size_cls0)) %>%
                mutate(
                    name_cls0 = factor(
                        name_cls0,
                        levels = unique(name_cls0),
                        ordered = TRUE
                    )
                )
            
        }
        
    )
    
)


IntercellCategoriesPairwise <- R6::R6Class(
    
    'IntercellCategoriesPairwise',
    
    inherit = Reader,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(fname = NULL){
            
            self$fname <- `if`(
                is.null(fname),
                omnipath2_settings$get(input_intercell_cat_pairwise),
                fname
            )
            
            super$initialize(name = self$fname)
            
            invisible(self)
            
        }
        
    )
    
)