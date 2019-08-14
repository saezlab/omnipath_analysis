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

require(tidyr)
require(dplyr)
require(rlang)
require(tibble)
require(R6)


ClusteringBase <- R6::R6Class(
    
    'ClusteringBase',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            data,
            key_var = NULL,
            row_var = NULL,
            value_var = NULL,
            hclust_args = list(),
            from_tidy = TRUE,
            ...
        ){
            
            self$data <- data
            self$hclust_args <- modifyList(
                list(
                    method = 'ward.D2'
                ),
                hclust_args
            )
            self$from_tidy <- from_tidy
            self$key_var <- enquo(key_var)
            self$value_var <- enquo(value_var)
            self$row_var <- enquo(row_var)
            
            self$main()
            
            invisible(self)
            
        },
        
        
        main = function(){
            
            private$preprocess()
            self$cluster()
            self$rearrange()
            
            invisible(self)
            
        },
        
        
        cluster = function(){
            
            self$clustering_x <- private$clustering(self$data_wide)
            self$clustering_y <- private$clustering(t(self$data_wide))
            self$order_x <- self$clustering_x$labels[self$clustering_x$order]
            self$order_y <- self$clustering_y$labels[self$clustering_y$order]
            
            invisible(self)
            
        },
        
        
        rearrange = function(){
            
            self$data <- self$get_ordered()
            
            invisible(self)
            
        },
        
        
        get_ordered = function(data = NULL){
            
            data <- `if`(is.null(data), self$data, data)
            
            (
                data %>%
                mutate(
                    !!self$key_var := factor(
                        !!self$key_var,
                        levels = self$order_y,
                        ordered = TRUE
                    ),
                    !!self$row_var := factor(
                        !!self$row_var,
                        levels = self$order_x,
                        ordered = TRUE
                    )
                ) %>%
                arrange(!!self$key_var, !!self$row_var)
            )
        
        }
        
        
    ),
    
    
    private = list(
        
        preprocess = function(){
            
            
            self$data_wide <- self$data %>%
                group_by(!!self$row_var, !!self$key_var) %>%
                summarize_all(first) %>%
                ungroup() %>%
                {`if`(
                    self$from_tidy,
                    select(
                        .,
                        !!self$key_var,
                        !!self$value_var,
                        !!self$row_var
                    ) %>%
                    spread(
                        key = !!self$key_var,
                        value = !!self$value_var
                    ) %>%
                    as.data.frame() %>%
                    column_to_rownames(quo_name(self$row_var)),
                    .
                )}
            
        },
        
        
        clustering = function(mat){
            
            do.call(
                hclust,
                c(
                    list(dist(mat)),
                    self$hclust_args
                )
            )
            
        }
        
    )
    
)
