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

require(dplyr)
require(ggplot2)
require(stringr)
require(R6)


IntercellClassSizes <- R6::R6Class(
    
    'IntercellClassSizes',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
                data,
                name,
                theme_args = list(),
                ...
            ){
            
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
                name = UQ(name),
                theme_args = theme_args,
                ...
            )
            
        },
        
        
        preprocess = function(...){
            
            self$data <- self$data %>%
                group_by(name0) %>%
                summarize_all(first) %>%
                mutate(
                    label0 = ifelse(is.na(label0), 'OmniPath', label0)
                ) %>%
                arrange(desc(label0 == 'OmniPath'), desc(size0)) %>%
                mutate(
                    label0 = factor(
                        label0, levels = unique(label0), ordered = TRUE
                    )
                )
            
            invisible(self)
            
        },
        
        
        plot = function(...){
            
            self$plt <- ggplot(self$data, aes(x = label0, y = size0)) +
                geom_col(fill = 'black') +
                xlab('Resources') +
                ylab('Number of proteins') +
                ggtitle(
                    gsub(
                        '_',
                        ' ',
                        first(str_to_title(self$data$parent0))
                    )
                )
            
            invisible(self)
            
        }
        
    )
    
)


IntercellClassSizesDots <- R6::R6Class(
    
    'IntercellClassSizesDots',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            data = NULL,
            ...
        ){
            
            self$input_param <- c('omnipath', 'proteins', 'all-class-levels')
            
            super$initialize(
                data = self$data,
                name = fig_cat_sizes_dots,
                width = 3.5,
                height = 4.5,
                theme_args = x_vertical_labels()
            )
            
            invisible(self)
            
        },
        
        
        plot = function(...){
            
            self$plt <- ggplot(
                    self$data,
                    aes(
                        x = cls,
                        y = database,
                        size = size,
                        color = database == 'OmniPath',
                        alpha = total_unique
                    )
                ) +
                geom_point(stroke = 0) +
                scale_color_manual(
                    values = omnipath2_settings$get(palette),
                    guide = FALSE
                ) +
                scale_alpha_manual(
                    values = c(
                        'total' = .66,
                        'unique' = 1.
                    ),
                    guide = FALSE
                ) +
                xlab('Inter-cellular\ncommunication roles') +
                ylab('Resources') +
                scale_size_area(
                    guide = guide_legend(
                        title = 'Number of\nproteins',
                        override.aes = list(
                            color = omnipath2_settings$get(palette)[1]
                        )
                    )
                )
            
            invisible(self)
            
        },
        
        
        preprocess = function(...){
            
            databases <- unique(self$data$resource0)
            
            classes <- unique(self$data$label0)

            self$by_entity <- self$by_entity %>%
                filter(
                    !is_complex &
                    !is.na(label)
                ) %>%
                group_by(entity_id, label) %>%
                mutate(is_unique = n() == 1) %>%
                ungroup() %>%
                filter(is_unique) %>%
                group_by(label, resource) %>%
                mutate(n_unique = n()) %>%
                summarize_all(first) %>%
                ungroup() %>%
                select(
                    database = resource,
                    cls = label,
                    size = n_unique
                ) %>%
                bind_rows(
                    {.} %>%
                    group_by(cls) %>%
                    filter(database != 'ComPPI' | cls != 'Extracellular') %>%
                    mutate(
                        size = sum(size),
                        database = 'OmniPath'
                    ) %>%
                    summarize_all(first) %>%
                    ungroup()
                )
            
            data_full <- expand.grid(
                    database = databases,
                    cls = classes
                ) %>%
                mutate(
                    database = as.character(database),
                    cls = as.character(cls)
                ) %>%
                left_join(
                    self$data %>%
                    select(label0, resource0, size0) %>%
                    group_by(label0, resource0) %>%
                    summarize_all(first) %>%
                    ungroup(),
                    by = c(
                        'database' = 'resource0',
                        'cls' = 'label0'
                    )
                ) %>%
                rename(size = size0) %>%
                mutate(
                    database = replace_na(database, 'OmniPath'),
                    total_unique = 'total'
                ) %>%
                filter(!is.na(size)) %>%
                bind_rows(
                    self$by_entity %>%
                    mutate(total_unique = 'unique')
                ) %>%
                arrange(desc(database == 'OmniPath'), database) %>%
                mutate(
                    database = factor(
                        database,
                        levels = unique(database),
                        ordered = TRUE
                    )
                )
            
            self$data <- data_full
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            self$data <- `if`(
                is.null(self$data),
                IntercellCategoriesPairwise$new(
                    input_param = self$input_param
                )$data,
                self$data
            )
            
            self$by_entity <- IntercellAnnotationByEntity$new()$data
            
            invisible(self)
            
        }
        
    )
    
)


IntercellClassSizesSeries <- R6::R6Class(
    
    'IntercellClassSizesSeries',
    
    inherit = CategoriesPairwisePlotSeries,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
                data = NULL,
                input_param = NULL,
                ...
            ){
            
            super$initialize(
                data = self$data,
                slice_var = parent0,
                plotter = IntercellClassSizes,
                name = quote(fig_cat_sizes),
                input_param = input_param,
                width_by = label0,
                width_min = .7,
                width_step = .23,
                ...
            )
            
            invisible(self)
            
        },
        
        
        preprocess = function(...){
            
            self$data <- self$data %>%
                filter(
                    name0 == name1 &
                    aspect0 == 'functional' &
                    source0 == 'composite' &
                    !is.na(name0)
                ) %>%
                {`if`(
                    'entity' %in% names(.),
                    filter(., entity == 'protein'),
                    .
                )}
            
            invisible(self)
            
        }
        
    )
    
)