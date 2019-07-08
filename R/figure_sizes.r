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
require(ggplot2)
require(R6)


CategorySizes <- R6::R6Class(
    
    'CategorySizes',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(data, name = NULL, ...){
            
            name <- `if`(is.null(name), 'sizes', name)
            
            super$initialize(
                data = data,
                name = name,
                ...
            )
            
        },
        
        
        preprocess = function(...){
            
            self$data <- self$data %>%
                group_by(name_cls0) %>%
                summarize_all(first) %>%
                mutate(
                    label0 = ifelse(is.na(label0), 'Total', label0)
                ) %>%
                arrange(desc(size_cls0)) %>%
                mutate(
                    label0 = factor(
                        label0, levels = unique(label0), ordered = TRUE
                    )
                )
            
            invisible(self)
            
        },
        
        
        plot = function(...){
            
            self$plt <- ggplot(self$data, aes(x = label0, y = size_cls0)) +
                geom_col(fill = 'black') +
                xlab('Resources') +
                ylab('Number of proteins') +
                theme(
                    axis.text.x = element_text(
                        angle = 90, vjust = 0.5, hjust = 1
                    )
                ) +
                ggtitle(
                    sprintf(
                        'Resources: %s',
                        toupper(first(self$data$parent0))
                    )
                )
            
            invisible(self)
            
        }
        
    )
    
)


CategorySizesSeries <- R6::R6Class(
    
    'CategorySizesSeries',
    
    inherit = CategoriesPairwisePlotSeries,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(data = NULL, ...){
            
            super$initialize(
                data = self$data,
                slice_var = parent0,
                plotter = CategorySizes,
                name = 'sizes',
                ...
            )
            
            invisible(self)
            
        },
        
        
        preprocess = function(...){
            
            self$data <- self$data %>%
                filter(
                    typ_cls0 == typ_cls1 &
                    entity == 'protein' &
                    typ_cls0 != 'misc' &
                    typ_cls0 != 'small_main' &
                    !is.na(name_cls0)
                )
            
            invisible(self)
            
        }
        
    )
    
)