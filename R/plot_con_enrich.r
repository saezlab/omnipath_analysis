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
            ...
        ){
            
            self$entity_type <- entity_type
            self$exclude <- exclude
            private$ensure_data(data)
            
            theme_args <- modifyList(
                list(
                    axis.text.x = element_text(
                        angle = 90, vjust = .5, hjust = 1, color = '#000000'
                    ),
                    panel.grid.major.y = element_blank()
                ),
                theme_args
            )
            
            super$initialize(
                data = self$data,
                name = 'connection-enrichment',
                theme_args = theme_args,
                width = 5,
                height = 4,
                ...
            )
            
        },
        
        plot = function(){
            
            self$plt <- ggplot(
                self$data,
                aes(x = cls_label0, y = cls_label1, fill = enrichment_log2)
            ) +
            geom_tile() +
            scale_fill_viridis_c(
                na.value = '#FFFFFF',
                guide = guide_colorbar(
                    title = 'Enrichment of\nconnections (log2)'
                ),
                limits = c(-2, 2)
            ) +
            ggtitle('Enrichment of network connections') +
            xlab('Inter-cellular\ncommunication roles') +
            ylab('Inter-cellular\ncommunication roles')
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            self$data <- self$data %>%
                filter(
                    !(typ_cls0 %in% self$exclude) &
                    !(typ_cls1 %in% self$exclude)
                ) %>%
                mutate(
                    con_expected =(
                        as.numeric(deg_total0) *
                        as.numeric(deg_total1) /
                        con_omnipath
                    ),
                    enrichment = con_all / con_expected
                ) %>%
                mutate(
                    enrichment_fold = ifelse(
                        enrichment > 1,
                        enrichment,
                        -1 / enrichment
                    )
                ) %>%
                mutate(
                    enrichment_log2 = ifelse(
                        abs(log2(enrichment)) > 2,
                        sign(enrichment) * 4,
                        enrichment
                    )
                ) %>%
                filter(!is.infinite(enrichment_log2))
            
            self$data <- bind_rows(
                self$data,
                self$data %>%
                    rename(
                        cls_label0 = cls_label1,
                        cls_label1 = cls_label0
                    ) %>%
                    filter(cls_label0 != cls_label1)
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
            
        }
        
    )
    
)
