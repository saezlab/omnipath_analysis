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
require(dplyr)
require(R6)


ConnectionEnrichment <- R6::R6Class(
    
    'ConnectionEnrichment',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            data = NULL,
            theme_args = list(),
            exclude = c('sub', 'misc'),
            entity_type = 'protein',
            cluster = TRUE,
            directed = FALSE,
            sign = NULL,
            order_by_undirected = TRUE,
            ...
        ){
            
            self$entity_type <- entity_type
            self$exclude <- exclude
            self$cluster <- cluster
            self$directed <- directed
            self$sign <- sign
            self$order_by_undirected <- order_by_undirected
            private$ensure_data(data)
            private$set_name()
            
            theme_args <- modifyList(
                c(
                    x_vertical_labels(),
                    list(
                        panel.grid.major.y = element_blank()
                    )
                ),
                theme_args
            )
            
            #return(private$setup())
            super$initialize(
                data = self$data,
                name = self$name,
                theme_args = theme_args,
                width = 5.075,
                height = 4,
                ...
            )
            
        },
        
        plot = function(){
            
            order_x <- unique(self$data$cls_label0)
            order_y <- unique(self$data$cls_label1)
            
            data <- self$data %>%
                arrange(cls_label0, cls_label1) %>%
                mutate(
                    cls_label0 = factor(
                        cls_label0,
                        levels = unique(cls_label0),
                        ordered = TRUE
                    ),
                    cls_label1 = factor(
                        cls_label1,
                        levels = unique(cls_label1),
                        ordered = TRUE
                    ),
                )
            
            data <- bind_rows(
                    private$upper_triangle(
                        data %>% filter(enrich_dir)
                    ),
                    private$upper_triangle(
                        data %>% filter(!enrich_dir),
                        swap_labels = TRUE
                    )
                ) %>%
                mutate(
                    sym = ifelse(enrich_dir, "\u25E4", "\u25E2")
                )
            
            self$plt <- ggplot(
                data,
                aes(
                    x = cls_label0,
                    y = cls_label1,
                    color = enrichment_log2,
                    label = sym
                )
            ) +
            geom_text(size = 4.2) +
            scale_color_viridis_c(
                na.value = '#FFFFFF',
                guide = guide_colorbar(
                    title = 'Enrichment of\nconnections (log2)'
                ),
                limits = c(-2, 2)
            ) +
            ggtitle(self$title) +
            xlab('Inter-cellular\ncommunication roles') +
            ylab('Inter-cellular\ncommunication roles')
            
        }
        
    ),
    #sym = ifelse(grp == 'a', "\u25E4", "\u25E2")
    
    private = list(
        
        setup = function(){
            
            private$add_enrichment()
            
            #return(private$set_order())
            private$set_order()
            
        },
        
        
        set_name = function(){
            
            self$name <- sprintf(
                'connection-enrichment-%s',
                private$dir_sign_choice(c('undir', 'dir', 'stim', 'inh'))
            )
            
            self$title <- sprintf(
                'Enrichment of network connections (%s)',
                private$dir_sign_choice(
                    c('all', 'directed', 'stimulation', 'inhibition')
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
            
            data <- private$get_enrichment()$data
            
            bind_rows(
                data %>% mutate(enrich_dir = TRUE),
                data %>% mutate(enrich_dir = FALSE)
            )
            
        },
        
        
        get_enrichment_directed = function(){
            
            bind_rows(
                private$get_enrichment(direction = c(0, 1))$data %>%
                    mutate(enrich_dir = TRUE),
                private$get_enrichment(direction = c(1, 0))$data %>%
                    mutate(enrich_dir = FALSE)
            )
            
        },
        
        
        set_variables = function(direction = c(0, 1)){
            
            dir_str <- do.call(
                sprintf,
                as.list(c('%sto%s', direction))
            )
            i <- private$dir_sign_choice(c(1, 2, 3, 4))
            
            self$variables <- list(
                side0 = list(
                    sym('deg_total0'),
                    sym(sprintf('deg_out%i', direction[1])),
                    sym(sprintf('deg_out%i_stim', direction[1])),
                    sym(sprintf('deg_out%i_inh', direction[1]))
                )[[i]],
                side1 = list(
                    sym('deg_total1'),
                    sym(sprintf('deg_in%i', direction[2])),
                    sym(sprintf('deg_in%i_stim', direction[2])),
                    sym(sprintf('deg_in%i_inh', direction[2]))
                )[[i]],
                observed = list(
                    sym('con_all'),
                    sym(sprintf('con_%s', dir_str)),
                    sym(sprintf('con_%s_stim', dir_str)),
                    sym(sprintf('con_%s_inh', dir_str))
                )[[i]],
                universe = list(
                    sym('con_omnipath'),
                    sym('con_omnipath_dir'),
                    sym('con_omnipath_stim'),
                    sym('con_omnipath_inh')
                )[[i]]
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
                    entity_type = self$entity_type
                )$data,
                data
            )
            
        },
        
        
        upper_triangle = function(data, swap_labels = FALSE){
            
            order_x <- unique(data$cls_label0)
            order_y <- unique(data$cls_label1)
            
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
                    match(cls_label0, order_x) <=
                    match(cls_label1, order_y)
                )
            )
            
        }
        
        
    )
    
)
