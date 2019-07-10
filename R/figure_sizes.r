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
        
        initialize = function(data, name = NULL, theme_args = list(), ...){
            
            name <- `if`(is.null(name), 'sizes', name)
            theme_args <- modifyList(
                list(
                    axis.text.x = element_text(
                        angle = 45, vjust = 1, hjust = 1, color = '#000000'
                    )
                ),
                theme_args
            )
            
            #private$set_width()
            
            super$initialize(
                data = data,
                name = name,
                theme_args = theme_args,
                ...
            )
            
        },
        
        
        preprocess = function(...){
            
            self$data <- self$data %>%
                group_by(name_cls0) %>%
                summarize_all(first) %>%
                mutate(
                    label0 = ifelse(is.na(src_label0), 'Total', src_label0)
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
