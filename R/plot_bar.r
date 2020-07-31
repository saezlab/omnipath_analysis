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

require(ggplot2)
require(ggnewscale)
require(dplyr)
require(R6)
require(scales)


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
            input_param,
            order_from_heatmap = TRUE,
            only_main_classes = TRUE,
            vertical = TRUE
        ){
            
            self$yvar <- enquo(yvar)
            self$order_from_heatmap <- order_from_heatmap
            self$ylab <- ylab
            self$title <- title
            self$vertical <- vertical
            self$input_param <- input_param
            self$only_main_classes <- only_main_classes
            
            size <- `if`(only_main_classes, 2.7, 4)
            
            super$initialize(
                data = NULL,
                name = UQ(enquo(name)),
                fname_param = c(
                    self$input_param,
                    list(
                        `if`(
                            self$only_main_classes,
                            'all-categories',
                            'main-categories'
                        )
                    )
                ),
                width = `if`(vertical, 2.4, size),
                height = `if`(vertical, size, 2.7),
                theme_args = `if`(vertical, NULL, x_vertical_labels())
            )
            
            invisible(self)
            
        },
        
        
        plot = function(...){
            
            self$plt <- ggplot(
                self$data,
                aes(x = label0, y = !!self$yvar)
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
                ConnectionEnrichment$new(
                    ordering_only = TRUE,
                    input_param = self$input_param,
                    only_main_classes = self$only_main_classes
                )$data,
                IntercellCategoriesPairwise$new(
                    input_param = self$input_param
                )$data
            ) %>%
            group_by(label0) %>%
            summarize_all(first) %>%
            ungroup()
            
            invisible(self)
            
        }
        
    )
    
)


StackedGroupedBarDot <- R6::R6Class(
    
    'StackedGroupedBarDot',
    
    lock_objects = FALSE,
    
    public = list(
        initialize = function(
            obj,
            data,
            xvar,
            yvar,
            fillvar,
            xlab,
            ylab,
            alphavar = NULL,
            log_y = FALSE,
            bar = TRUE,
            color_values = NULL,
            color_labels = NULL,
            legend_title = NULL,
            shape = 16,
            size = 5,
            scale_alpha_param = NULL,
            position = 'identity',
            expand_y = waiver(),
            tick_label_size = 8
        ){
            
            self$obj <- obj
            self$data <- data
            self$xvar <- enquo(xvar)
            self$yvar <- enquo(yvar)
            self$fillvar <- enquo(fillvar)
            self$alphavar <- enquo(alphavar)
            self$xlab <- xlab
            self$ylab <- ylab
            self$log_y <- log_y
            self$bar <- bar
            self$color_values <- color_values
            self$color_labels <- color_labels
            self$legend_title <- legend_title
            self$shape <- shape
            self$size <- size
            self$scale_alpha_param <- scale_alpha_param
            self$position <- position
            self$expand_y <- expand_y
            self$tick_label_size <- tick_label_size
            
            self$main()
            
            invisible(self)
            
        },
        
        main = function(){
            
            self$obj$plt <- ggplot(
                    self$data,
                    aes(
                        y = !!self$yvar,
                        x = !!self$xvar
                    )
                ) %>%
                `+`(
                    {`if`(
                        self$bar,
                        geom_col(
                            aes(fill = !!self$fillvar),
                            position = `if`(self$log_y, 'dodge', 'stack')
                        ),
                        geom_point(
                            aes(
                                color = !!self$fillvar,
                                alpha = !!self$alphavar
                            ),
                            size = self$size,
                            shape = self$shape,
                            position = self$position
                        )
                    )}
                ) %>%
                `+`(
                    {`if`(
                        self$bar,
                        scale_fill_manual,
                        scale_color_manual
                    )(
                        values = self$color_values,
                        labels = self$color_labels,
                        name = self$legend_title
                    )}
                ) %>%
                `+`(
                    {`if`(
                        quo_text(self$alphavar) == 'NULL',
                        NULL,
                        do.call(
                            scale_alpha_manual,
                            self$scale_alpha_param
                        )
                    )}
                ) %>%
                `+`(
                    {`if`(
                        self$log_y,
                        scale_y_log10(labels = comma, expand = self$expand_y),
                        NULL
                    )}
                ) %>%
                `+`(coord_flip()) %>%
                `+`(xlab(self$xlab)) %>%
                `+`(ylab(self$ylab)) %>%
                `+`(
                    theme(
                        axis.text.x = element_text(size = self$tick_label_size),
                        axis.text.y = element_text(size = self$tick_label_size)
                    )
                )
            
            invisible(self)
            
        }
        
    )
    
)


IntercellSizesBar <- R6::R6Class(
    
    'IntercellSizesBar',
    
    inherit = SimpleBar,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            input_param,
            only_main_classes = TRUE,
            ...
        ){
            
            super$initialize(
                yvar = size_log10,
                name = fig_intercell_size,
                ylab = 'Number of proteins',
                title = paste0(
                    'Number of proteins\n',
                    'in inter-cellular communication roles'
                ),
                order_from_heatmap = TRUE,
                vertical = TRUE,
                input_param = input_param,
                only_main_classes = only_main_classes
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(...){
            
            super$setup()
            
            self$data <- self$data %>%
                mutate(size_log10 = log10(size0))
            
            invisible(self)
            
        }
        
    )
    
)


IntercellCoverageBar <- R6::R6Class(
    
    'IntercellCoverageBar',
    
    inherit = SimpleBar,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            input_param,
            only_main_classes = TRUE,
            ...
        ){
            
            super$initialize(
                yvar = omnipath_coverage,
                name = fig_intercell_cov,
                ylab = 'Coverage',
                title = paste0(
                    'Coverage by OmniPath:\n',
                    'inter-cellular communication roles'
                ),
                order_from_heatmap = TRUE,
                vertical = FALSE,
                input_param = input_param,
                only_main_classes = only_main_classes
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(...){
            
            super$setup()
            
            self$data <- self$data %>%
                mutate(omnipath_coverage = in_network_cls0 / size0 * 100)
            
            invisible(self)
            
        }
        
    )
    
)


IntercellDegreeBar <- R6::R6Class(
    
    'IntercellDegreeBar',
    
    inherit = SimpleBar,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            input_param,
            only_main_classes = TRUE,
            ...
        ){
            
            super$initialize(
                yvar = mean_degree,
                name = fig_intercell_deg,
                ylab = 'Degree',
                title = paste0(
                    'Mean degree in OmniPath:\n',
                    'inter-cellular communication roles'
                ),
                order_from_heatmap = TRUE,
                vertical = TRUE,
                input_param = input_param,
                only_main_classes = only_main_classes
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(...){
            
            super$setup()
            
            self$data <- self$data %>%
                mutate(mean_degree = deg_total0 / size0)
            
            invisible(self)
            
        }
        
    )
    
)
