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


Enrichment <- R6::R6Class(
    
    'Enrichment',
    
    lock_objects = FALSE,
    
    public = list(
        
        
        initialize = function(
            data,
            side0,
            side1,
            observed,
            universe,
            directed = TRUE,
            exclude = NULL,
            enrich_var = NULL,
            expected_var = NULL
        ){
            
            self$data <- data
            self$side0 <- enquo(side0)
            self$side1 <- enquo(side1)
            self$observed <- enquo(observed)
            self$universe <- enquo(universe)
            self$directed <- directed
            self$exclude <- exclude
            self$enrich_var <- `if`(
                is.null(enrich_var),
                sym('enrichment'),
                enquo(enrich_var)
            )
            self$enrich_fold_var <- sym(
                sprintf('%s_fold', quo_text(self$enrich_var))
            )
            self$enrich_log2_var <- sym(
                sprintf('%s_log2', quo_text(self$enrich_var))
            )
            self$expected_var <- `if`(
                is.null(expected_var),
                sym('expected'),
                enquo(expected_var)
            )
            
            self$main()
            
            invisible(self)
            
        },
        
        
        main = function(){
            
            self$filter()
            self$enrichment()
            self$add_reverse_direction()
            
            invisible(self)
            
        },
        
        
        filter = function(){
            
            self$data <- self$data %>%
                filter(
                    !(typ_cls0 %in% self$exclude) &
                    !(typ_cls1 %in% self$exclude)
                )
            
            invisible(self)
            
        },
        
        
        enrichment = function(){
            
            self$data <- self$data %>%
                mutate(
                    !!self$expected_var := (
                        as.numeric(!!self$side0) *
                        as.numeric(!!self$side1) /
                        !!self$observed
                    )
                )%>%
                mutate(
                    !!self$enrich_var := (
                        !!self$universe /
                        !!self$expected_var
                    )
                ) %>%
                mutate(
                    !!self$enrich_fold_var := ifelse(
                        !!self$enrich_var > 1,
                        !!self$enrich_var,
                        -1 / !!self$enrich_var
                    ),
                    !!self$enrich_log2_var := log2(!!self$enrich_var)
                ) %>%
                mutate(
                    !!self$enrich_log2_var := ifelse(
                        abs(!!self$enrich_log2_var) > 2,
                        sign(!!self$enrich_log2_var) * 2,
                        !!self$enrich_log2_var
                    )
                ) %>%
                filter(!is.infinite(enrichment_log2))
            
            invisible(self)
            
        },
        
        
        add_reverse_direction = function(){
            
            if(!self$directed){
                
                self$data <- bind_rows(
                    self$data,
                    self$data %>%
                        rename(
                            cls_label0 = cls_label1,
                            cls_label1 = cls_label0
                        ) %>%
                        filter(cls_label0 != cls_label1)
                )
                
            }
            
            invisible(self)
            
        }
        
        
    )
    
)
