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

require(ggplot2)
require(ggnewscale)
require(dplyr)
require(R6)


SimpleBar <- R6::R6Class(
    
    'SimpleBar',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            yvar,
            name,
            ylab,
            title,
            order_from_heatmap = TRUE,
            vertical = TRUE
        ){
            
            self$yvar <- enquo(yvar)
            self$order_from_heatmap <- order_from_heatmap
            self$ylab <- ylab
            self$title <- title
            self$vertical <- vertical
            
            super$initialize(
                data = NULL,
                name = name,
                width = `if`(vertical, 2.4, 4),
                height = `if`(vertical, 4, 2.7),
                theme_args = `if`(vertical, NULL, x_vertical_labels())
            )
            
            invisible(self)
            
        },
        
        
        plot = function(...){
            
            self$plt <- ggplot(
                self$data,
                aes(x = cls_label0, y = !!self$yvar)
            ) %>%
            `+`(geom_col(fill = '#333333')) %>%
            {`if`(
                self$vertical,
                `+`(`+`(., coord_flip()), scale_y_reverse()),
                .
            )} %>%
            `+`(ylab(self$ylab)) %>%
            `+`(xlab('Inter-cellular\ncommunication roles')) %>%
            `+`(ggtitle(self$title))
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            self$data <- `if`(
                self$order_from_heatmap,
                ConnectionEnrichment$new(ordering_only = TRUE)$data,
                IntercellCategoriesPairwise$new()$data
            ) %>%
            group_by(cls_label0) %>%
            summarize_all(first) %>%
            ungroup()
            
            invisible(self)
            
        }
        
    )
    
)


SizesBar <- R6::R6Class(
    
    'SizesBar',
    
    inherit = SimpleBar,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                yvar = size_log10,
                name = 'intercell-sizes',
                ylab = 'Number of proteins',
                title = paste0(
                    'Number of proteins\n',
                    'in inter-cellular communication roles'
                ),
                order_from_heatmap = TRUE,
                vertical = TRUE
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(...){
            
            super$setup()
            
            self$data <- self$data %>%
                mutate(size_log10 = log10(size_cls0))
            
            invisible(self)
            
        }
        
    )
    
)


CoverageBar <- R6::R6Class(
    
    'CoverageBar',
    
    inherit = SimpleBar,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                yvar = omnipath_coverage,
                name = 'intercell-cov',
                ylab = 'Coverage',
                title = paste0(
                    'Coverage by OmniPath:\n',
                    'inter-cellular communication roles'
                ),
                order_from_heatmap = TRUE,
                vertical = FALSE
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(...){
            
            super$setup()
            
            self$data <- self$data %>%
                mutate(omnipath_coverage = in_omnipath_cls0 / size_cls0 * 100)
            
            invisible(self)
            
        }
        
    )
    
)


DegreeBar <- R6::R6Class(
    
    'DegreeBar',
    
    inherit = SimpleBar,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                yvar = mean_degree,
                name = 'intercell-deg',
                ylab = 'Degree',
                title = paste0(
                    'Mean degree in OmniPath:\n',
                    'inter-cellular communication roles'
                ),
                order_from_heatmap = TRUE,
                vertical = TRUE
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(...){
            
            super$setup()
            
            self$data <- self$data %>%
                mutate(mean_degree = deg_total0 / size_cls0)
            
            invisible(self)
            
        }
        
    )
    
)
