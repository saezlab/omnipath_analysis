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
require(R6)


ResourceCoverage <- R6::R6Class(
    
    'ResourceCoverage',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            con_enrich_input_param,
            res_by_entity_input_param,
            ann_by_entity_input_param = list(),
            only_main_classes = FALSE
        ){
            
            self$con_enrich_input_param <- con_enrich_input_param
            self$res_by_entity_input_param <- res_by_entity_input_param
            self$ann_by_entity_input_param <- ann_by_entity_input_param
            self$only_main_classes <- only_main_classes
            
            super$initialize(
                data = self$data,
                name = fig_res_cov,
                fname_param = c(
                    self$con_enrich_input_param,
                    self$res_by_entity_input_param,
                    self$ann_by_entity_input_param,
                    list(
                        `if`(
                            only_main_classes,
                            'main-classes',
                            'all-classes'
                        )
                    )
                ),
                height = 7,
                width = 4.5,
                theme_args = c(
                    x_vertical_labels(),
                    list(
                        panel.grid = element_blank()
                    )
                )
            )
            
            invisible(self)
            
        },
        
        
        plot = function(...){
            
            self$plt <- ggplot(
                self$data,
                aes(x = class_label, y = resource, fill = log10(coverage))
            ) +
            geom_tile() +
            scale_fill_viridis_c(
                guide = guide_colorbar(
                    title = 'Coverage\n[log10(%)]'
                )
            ) +
            xlab('Inter-cellular\ncommunication roles') +
            ylab('Network resources') +
            ggtitle(paste0(
                'Coverage of network resources\n',
                'on inter-cellular signaling proteins'
            ))
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            self$order <- (
                ConnectionEnrichment$new(
                    ordering_only = TRUE,
                    input_param = self$con_enrich_input_param,
                    only_main_classes = self$only_main_classes
                )$clustering$order_x
            )
            
            res <- ResourceByEntity$new(
                    input_param = self$res_by_entity_input_param
                )$data %>%
                filter(!is_complex) %>%
                select(entity_id, resource)
            
            icc <- IntercellAnnotationByEntity$new(
                    input_param = self$ann_by_entity_input_param
                )$data %>%
                {`if`(
                    self$only_main_classes,
                    filter(
                        .,
                        cls %in% omnipath2_settings$get(
                            intercell_main_classes
                        )
                    ),
                    .
                )} %>%
                filter(!is_complex) %>%
                select(entity_id, class_label) %>%
                group_by(class_label) %>%
                mutate(class_size = n()) %>%
                ungroup() %>%
                inner_join(res, by = c('entity_id')) %>%
                group_by(resource, class_label) %>%
                mutate(in_resource = n()) %>%
                summarize_all(first) %>%
                select(-entity_id) %>%
                mutate(coverage = in_resource / class_size * 100) %>%
                ungroup() %>%
                filter(class_label %in% self$order) %>%
                mutate(
                    class_label = factor(
                        class_label,
                        levels = self$order,
                        ordered = TRUE
                    )
                )
            
            self$data <- icc
            
            self$height <- length(unique(icc$resource)) / 5
            self$width <- length(unique(icc$class_label)) / 6.6 + 2.6
            
            invisible(self)
            
        }
        
    )
    
)
