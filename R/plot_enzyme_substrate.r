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


EnzymeSubstrateBase <- R6::R6Class(
    
    'EnzymeSubstrateBase',
    
    lock_objects = FALSE,
    
    inherit = SinglePlot,
    
    public = list(
        
        initialize = function(
                name,
                complexes = FALSE,
                keep_references = FALSE
            ){
            
            self$complexes <- complexes
            self$keep_references <- keep_references
            theme_args <- x_vertical_labels()
            
            super$initialize(
                data = self$data,
                name = UQ(enquo(name)),
                theme_args = theme_args,
                width = 9
            )
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            self$enz_sub <- EnzymeSubstrate$new()$data
            
            if(self$complexes){
                
                complexes <- {
                    ComplexesByResource$new(expand_members = TRUE)$data
                } %>%
                mutate(complex_id = as.integer(as.factor(complex_id))) %>%
                select(-resource, -n_members)
                
            }
            
            self$enz_sub <- self$enz_sub %>%
                select(
                    -enzyme_genesymbol,
                    -substrate_genesymbol,
                    -isoforms
                ) %>%
                {`if`(
                    self$keep_references,
                    .,
                    select(., -references)
                )} %>%
                {`if`(
                    self$complexes,
                    left_join(., complexes, by = c('enzyme' = 'members')) %>%
                    left_join(
                        complexes,
                        by = c('substrate' = 'members'),
                        suffix = c('.enz', '.sub')
                    ) %>%
                    mutate(
                        complex_id.enz = replace_na(complex_id.enz, 0),
                        complex_id.sub = replace_na(complex_id.sub, 0),
                        in_complex = (
                            complex_id.enz != 0 &
                            complex_id.enz == complex_id.sub
                        )
                    ) %>%
                    group_by(
                        enzyme,
                        substrate,
                        residue_type,
                        residue_offset,
                        modification
                    ) %>%
                    mutate(in_complex = any(in_complex)) %>%
                    summarize_all(first) %>%
                    ungroup() %>%
                    select(
                        -complex_id.enz,
                        -complex_id.sub,
                    ) %>%
                    mutate(self_modification = enzyme == substrate) %>%
                    mutate(
                        category = ifelse(
                            self_modification,
                            'self',
                            ifelse(
                                in_complex,
                                'in_complex',
                                'between_entities'
                            )
                        )
                    ),
                    .
                )}
            
            
            self$enz_sub <- self$enz_sub %>%
                mutate(
                    shared = sapply(
                        str_split(sources, ';'),
                        function(s){
                            s1 <- sapply(
                                str_split(s, '_'),
                                function(s0){tail(s0, n = 1)}
                            )
                            length(unique(s1)) > 1
                        }
                    )
                )
            
            self$enz_sub_by_resource <- bind_rows(
                self$enz_sub %>%
                separate_rows(sources, sep = ';') %>%
                group_by(sources) %>%
                mutate(n_total = n()) %>%
                ungroup() %>%
                group_by(sources, shared) %>%
                mutate(n_shared = n()) %>%
                ungroup() %>%
                {`if`(
                    self$complexes,
                    group_by(., category, sources) %>%
                    mutate(n_by_category = n()) %>%
                    ungroup(),
                    .
                )},
                self$enz_sub %>%
                mutate(sources = 'Total') %>%
                mutate(n_total = n()) %>%
                group_by(shared) %>%
                mutate(n_shared = n()) %>%
                ungroup() %>%
                {`if`(
                    self$complexes,
                    group_by(., category) %>%
                    mutate(n_by_category = n()) %>%
                    ungroup(),
                    .
                )}
            ) %>%
            arrange(desc(n_total)) %>%
            mutate(
                sources = factor(
                    sources,
                    levels = unique(sources),
                    ordered = TRUE
                )
            )
            
            invisible(self)
            
        }
        
    )
    
)


EnzymeSubstrateShared <- R6::R6Class(
    
    'EnzymeSubstrateShared',
    
    inherit = EnzymeSubstrateBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(name = fig_enzyme_substrate_self)
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(
                    self$enz_sub_by_resource,
                    aes(x = sources, y = n_shared)
                ) +
                geom_col(aes(fill = shared)) +
                scale_fill_manual(
                    values = c(
                        `TRUE` = '#4268B3',
                        `FALSE` = '#B3C5E9'
                    ),
                    labels = c(
                        `TRUE` = 'Shared',
                        `FALSE` = 'Unique'
                    ),
                    name = 'Target of the\nenzyme-substrate\ninteraction'
                ) +
                xlab('Resources') +
                ylab('Enzyme-substrate\ninteractions')
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$enz_sub_by_resource <- self$enz_sub_by_resource %>%
                group_by(sources, shared) %>%
                summarize_all(first) %>%
                ungroup()
            
        }
        
    )
    
)


EnzymeSubstrateSelf <- R6::R6Class(
    
    'EnzymeSubstrateSelf',
    
    inherit = EnzymeSubstrateBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_enzyme_substrate_self,
                complexes = TRUE
            )
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(
                    self$enz_sub_by_resource,
                    aes(x = sources, y = n_by_category)
                ) +
                geom_col(aes(fill = category)) +
                scale_fill_manual(
                    values = c(
                        self = '#4268B3',
                        in_complex = '#6F8DCF',
                        between_entities = '#B3C5E9'
                    ),
                    labels = c(
                        self = 'Self',
                        in_complex = 'Within complex',
                        between_entities = 'Other entity'
                    ),
                    name = 'Enzyme-substrate\ninteractions'
                ) +
                xlab('Resources') +
                ylab('Enzyme-substrate\ninteractions')
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$enz_sub_by_resource <- self$enz_sub_by_resource %>%
                group_by(sources, category) %>%
                summarize_all(first) %>%
                ungroup()
            
        }
        
    )
    
)
