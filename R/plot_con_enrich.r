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


ConnectionEnrichment <- R6::R6Class(
    
    'ConnectionEnrichment',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            data = NULL,
            name = fig_con_enrich,
            input_param = NULL,
            fname_param = list(),
            theme_args = list(),
            exclude = c('sub', 'misc'),
            entity_type = 'protein',
            cluster = TRUE,
            directed = FALSE,
            sign = NULL,
            order_by_undirected = TRUE,
            ordering_only = FALSE,
            only_main_classes = FALSE,
            heatmap_variables = 'enrich-count',
            ...
        ){
            
            self$name <- enquo(name)
            self$entity_type <- entity_type
            self$exclude <- exclude
            self$cluster <- cluster
            self$directed <- directed
            self$sign <- sign
            self$order_by_undirected <- order_by_undirected
            self$ordering_only <- ordering_only
            self$only_main_classes <- only_main_classes
            self$heatmap_variables <- heatmap_variables
            self$input_param <- input_param
            private$ensure_data(data)
            private$set_name()
            
            theme_args <- modifyList(
                c(
                    x_vertical_labels(),
                    list(
                        panel.grid.major.y = element_blank(),
                        legend.title = element_text(size = 10),
                        legend.justification = 'bottom',
                        legend.text = element_text(size = 7, hjust = 1)
                    )
                ),
                theme_args
            )
            
            super$initialize(
                data = self$data,
                name = UQ(self$name),
                fname_param = self$fname_param,
                theme_args = theme_args,
                width = 5.075,
                height = 4,
                ...
            )
            
        },
        
        
        main = function(){
            
            if(self$ordering_only){
                
                private$setup()
                
            }else{
                
                super$main()
                
            }
            
        },
        
        
        plot = function(){
            
            data <- bind_rows(
                    private$upper_triangle(
                        self$data %>% filter(enrich_dir)
                    ),
                    private$upper_triangle(
                        self$data %>% filter(!enrich_dir),
                        swap_labels = TRUE
                    ) %>%
                    filter(cls_label0 != cls_label1)
                ) %>%
                mutate(
                    sym1 = ifelse(enrich_dir, "\u25E4", "\u25E2"),
                    sym2 = ifelse(enrich_dir, "\u25E2", "\u25E4")
                )
            
            self$plt <- ggplot(
                mapping = aes(
                    x = cls_label0,
                    y = cls_label1,
                    label = sym1
                )
            ) +
            geom_text(
                data = data,
                mapping = aes(
                    color = enrichment_log2
                ),
                size = 4.2
            ) +
            scale_color_viridis_c(
                na.value = '#FFFFFF',
                guide = guide_colorbar(
                    title = 'Enrichment of\nconnections (log2)',
                    barheight = unit(.5, 'inches')
                ),
                limits = c(-2, 2)
            ) +
            new_scale_color() +
            geom_text(
                data = data,
                mapping = aes(
                    x = cls_label1,
                    y = cls_label0,
                    color = numof_con_log10,
                    label = sym2
                )
            ) +
            scale_color_gradient(
                low = '#FFFFFF',
                high = '#333333',
                na.value = '#FFFFFF',
                limits = c(0, 5),
                guide = guide_colorbar(
                    title = 'Number of\nconnections (log10)',
                    barheight = unit(.5, 'inches')
                )
            ) +
            ggtitle(self$title) +
            xlab('Inter-cellular\ncommunication roles') +
            ylab('Inter-cellular\ncommunication roles')
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            private$add_enrichment()
            private$set_order()
            
        },
        
        
        set_name = function(){
            
            dir_sign <- private$dir_sign_choice(
                c('undir', 'dir', 'stim', 'inh')
            )
            
            self$fname_param <- c(
                self$input_param,
                list(
                    dir_sign,
                    `if`(
                        self$only_main_classes,
                        'main-categories',
                        'all-categories'
                    ),
                    self$heatmap_variables
                )
            )
            
            self$title <- sprintf(
                'Enrichment of network connections (%s)',
                dir_sign
            )
            
        },
        
        
        add_numof_con = function(data, direction = c(0, 1)){
            
            dir_str <- do.call(
                sprintf,
                as.list(c('%sto%s', direction))
            )
            
            var <- private$dir_sign_choice(c(
                sym('con_all'),
                sym(sprintf('con_%s', dir_str)),
                sym(sprintf('con_%s_stim', dir_str)),
                sym(sprintf('con_%s_inh', dir_str))
            ))
            
            (
                data %>%
                mutate(
                    numof_con = !!var,
                    numof_con_log10 = log10(!!var)
                )
            )
            
        },
        
        
        add_enrichment = function(){
            
            self$data <- `if`(
                self$directed,
                private$get_enrichment_directed(),
                private$get_enrichment_undirected()
            )
            
        },
        
        
        get_enrichment_undirected = function(){
            
            data <- private$get_enrichment()$data %>%
                private$add_numof_con()
            
            bind_rows(
                data %>% mutate(enrich_dir = TRUE),
                data %>% mutate(enrich_dir = FALSE)
            )
            
        },
        
        
        get_enrichment_directed = function(){
            
            data <- bind_rows(
                private$get_enrichment(direction = c(0, 1))$data %>%
                    private$add_numof_con(direction = c(0, 1)) %>%
                    mutate(enrich_dir = TRUE),
                private$get_enrichment(direction = c(1, 0))$data %>%
                    private$add_numof_con(direction = c(1, 0)) %>%
                    mutate(enrich_dir = FALSE)
            )
            
            bind_rows(
                data,
                data %>%
                    rename(
                        cls_label0 = cls_label1,
                        cls_label1 = cls_label0
                    ) %>%
                    filter(cls_label0 != cls_label1)
            )
            
        },
        
        
        set_variables = function(direction = c(0, 1)){
            
            dir_str <- do.call(
                sprintf,
                as.list(c('%sto%s', direction))
            )
            
            self$variables <- list(
                side0 = private$dir_sign_choice(c(
                    sym('deg_total0'),
                    sym(sprintf('deg_out%i', direction[1])),
                    sym(sprintf('deg_out%i_stim', direction[1])),
                    sym(sprintf('deg_out%i_inh', direction[1]))
                )),
                side1 = private$dir_sign_choice(c(
                    sym('deg_total1'),
                    sym(sprintf('deg_in%i', direction[2])),
                    sym(sprintf('deg_in%i_stim', direction[2])),
                    sym(sprintf('deg_in%i_inh', direction[2]))
                )),
                observed = private$dir_sign_choice(c(
                    sym('con_all'),
                    sym(sprintf('con_%s', dir_str)),
                    sym(sprintf('con_%s_stim', dir_str)),
                    sym(sprintf('con_%s_inh', dir_str))
                )),
                universe = private$dir_sign_choice(c(
                    sym('con_network'),
                    sym('con_network_dir'),
                    sym('con_network_stim'),
                    sym('con_network_inh')
                ))
            )
            
        },
        
        
        dir_sign_choice = function(choices){
            
            if(is.null(names(choices))){
                
                names(choices) <- c('undir', 'dir', 'stim', 'inh')
                
            }
            
            result <- `if`(
                self$directed,
                `if`(
                    is.null(self$sign),
                    choices[['dir']],
                    choices[[self$sign]]
                ),
                choices[['undir']]
            )
            
        },
        
        
        set_order = function(){
            
            if(self$cluster){
                
                directed <- self$directed
                self$directed <- !self$order_by_undirected && self$directed
                
                enrich <- private$get_enrichment()
                #return(enrich)
                self$clustering <- ClusteringBase$new(
                    enrich$data,
                    cls_label0,
                    cls_label1,
                    enrichment_log2
                )
                self$data <- self$clustering$get_ordered(data = self$data)
                
                self$directed <- directed
                
            }
            
        },
        
        
        get_enrichment = function(direction = c(0, 1)){
            
            private$set_variables(direction = direction)
            
            variables <- c(
                self$variables,
                list(
                    data = self$data,
                    directed = self$directed,
                    exclude = self$exclude
                )
            )
            
            do.call(
                Enrichment$new,
                variables
            )
            
        },
        
        
        ensure_data = function(data){
            
            self$data <- `if`(
                is.null(data),
                IntercellCategoriesPairwise$new(
                    input_param = self$input_param,
                    entity_type = self$entity_type
                )$data,
                data
            )
            private$filter_classes()
            
            invisible(self)
            
        },
        
        
        filter_classes = function(){
            
            main_classes <- omnipath2_settings$get(intercell_main_classes)
            
            self$data <- self$data %>%
                {`if`(
                    self$only_main_classes,
                    filter(
                        .,
                        cls0 %in% main_classes &
                        cls1 %in% main_classes
                    ),
                    .
                )}
            
            invisible(self)
            
        }
        
        
        upper_triangle = function(...){
            
            private$triangle(...)
            
        },
        
        
        lower_triangle = function(...){
            
            private$triangle(..., op = `>=`)
            
        },
        
        
        triangle = function(data, swap_labels = FALSE, op = `<=`){
            
            order_x <- levels(data$cls_label0)
            order_y <- levels(data$cls_label1)
            
            #return(data)
            
            (
                data %>%
                {`if`(
                    swap_labels,
                    rename(
                        .,
                        cls_label0 = cls_label1,
                        cls_label1 = cls_label0
                    ),
                    .
                )} %>%
                filter(
                    op(
                        match(cls_label0, order_x),
                        match(cls_label1, order_y)
                    )
                ) %>%
                mutate(
                    cls_label0 = factor(
                        cls_label0,
                        levels = order_x,
                        ordered = TRUE
                    ),
                    cls_label1 = factor(
                        cls_label1,
                        levels = order_y,
                        ordered = TRUE
                    )
                )
            )
            
        }
        
        
    )
    
)
