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
require(dplyr)
require(stringr)
require(R6)


EnzymeSubstrateBase <- R6::R6Class(
    
    'EnzymeSubstrateBase',
    
    lock_objects = FALSE,
    
    inherit = SinglePlot,
    
    public = list(
        
        initialize = function(
                name,
                complexes = FALSE,
                keep_references = FALSE,
                xlab_vertical = TRUE,
                width = 9,
                ...
            ){
            
            self$complexes <- complexes
            self$keep_references <- keep_references
            theme_args <- `if`(xlab_vertical, x_vertical_labels(), list())
            args <- modifyList(
                list(
                    data = self$data,
                    name = enquo(name),
                    theme_args = theme_args,
                    width = width
                ),
                list(...)
            )
            
            do.call(
                super$initialize,
                args
            )
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            self$enz_sub <- EnzymeSubstrate$new()$data
            
            if(self$complexes){
                
                complexes <- {
                    Complexes$new(
                        expand_members = TRUE,
                        keep_references = FALSE
                    )$data
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
                group_by(sources, modification) %>%
                mutate(n_modtype = n()) %>%
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
                mutate(sources = 'OmniPath') %>%
                mutate(n_total = n()) %>%
                group_by(shared) %>%
                mutate(n_shared = n()) %>%
                ungroup() %>%
                group_by(modification) %>%
                mutate(n_modtype = n()) %>%
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
            
            super$initialize(name = fig_enzyme_substrate_shared)
            
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
        
        initialize = function(name = NULL, log_y = FALSE, bar = TRUE){
            
            self$log_y <- log_y
            self$bar <- bar
            
            super$initialize(
                name = fig_enzyme_substrate_self,
                fname_param = list(
                    `if`(self$log_y, 'dodge', 'stack'),
                    `if`(self$bar, 'bar', 'dot')
                ),
                complexes = TRUE,
                width = 6,
                height = 6,
                xlab_vertical = FALSE
            )
            
            gc()
            
            invisible(self)
            
        },
        
        
        plot = function(){
            
            StackedGroupedBarDot$new(
                obj = self,
                data = self$enz_sub_by_resource,
                xvar = sources,
                yvar = n_by_category,
                fillvar = category,
                xlab = 'Resources',
                ylab = sprintf(
                    'Enzyme-substrate\ninteractions%s',
                    `if`(self$log_y, ' (log)', '')
                ),
                log_y = self$log_y,
                bar = self$bar,
                color_values = `names<-`(
                    omnipath2_settings$get(three_shades_1),
                    c('self', 'in_complex', 'between_entities')
                ),
                color_labels = c(
                        self = 'Self',
                        in_complex = 'Within complex',
                        between_entities = 'Other entity'
                    ),
                legend_title = 'Target of the\nenzyme-substrate\ninteractions'
            )
            
            invisible(self)
            
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


EnzymeSubstrateModtype <- R6::R6Class(
    
    'EnzymeSubstrateModtype',
    
    inherit = EnzymeSubstrateBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(name = NULL, ...){
            
            name <- `if`(
                'name' %in% names(match.call()),
                enquo(name),
                quo(fig_enzyme_substrate_modtype)
            )
            
            args <- modifyList(
                list(
                    name = name,
                    width = 9,
                    height = 7
                ),
                list(...)
            )
            
            do.call(super$initialize, args)
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(
                    self$enz_sub_by_resource,
                    aes(x = sources, y = n_modtype)
                ) +
                geom_col(aes(fill = modification)) +
                scale_fill_discrete(
                    guide = guide_legend(
                        title = 'Modification type'
                    )
                ) +
                xlab('Resources') +
                ylab('Enzyme-substrate\ninteractions')
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$enz_sub_by_resource <- self$enz_sub_by_resource %>%
                group_by(sources, modification) %>%
                summarize_all(first) %>%
                ungroup() %>%
                mutate(modification = str_to_title(modification)) %>%
                filter(!is.na(modification) & n_modtype > 10) %>%
                # it does not look good by this ordering
                arrange(desc(n_modtype)) %>%
                mutate(
                    modification = factor(
                        modification,
                        rev(unique(modification)),
                        ordered = TRUE
                    )
                )
            
            invisible(self)
            
        }
        
    )
    
)


EnzymeSubstrateModtypeDot <- R6::R6Class(
    
    'EnzymeSubstrateModtypeDot',
    
    inherit = EnzymeSubstrateModtype,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_enzyme_substrate_modtype_dot,
                height = 5,
                xlab_vertical = FALSE
            )
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(
                    self$enz_sub_by_resource,
                    aes(y = sources, x = n_modtype, color = modification)
                ) +
                geom_point(size = 5, alpha = .7) +
                scale_color_discrete(
                    guide = guide_legend(
                        title = 'Modification type'
                    )
                ) +
                scale_x_log10() +
                ylab('Resources') +
                xlab('Enzyme-substrate interactions')
            
            invisible(self)
            
        }
        
    )
    
)


EnzymeSubstrateNumofResources <- R6::R6Class(
    
    'EnzymeSubstrateNumofResources',
    
    inherit = EnzymeSubstrateBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_enzyme_substrate_numof_res,
                height = 4,
                width = 5,
                xlab_vertical = FALSE,
                theme_args = list(
                    axis.text.x = element_text(size = 14)
                )
            )
            
        },
        
        
        plot = function(){
            
            self$plt <- ggplot(
                    self$enz_sub,
                    aes(x = n_resources)
                ) +
                geom_bar(fill = 'black', stat = 'count') +
                xlab('Number of resources') +
                ylab('Enzyme-substrate interactions')
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$enz_sub <- self$enz_sub %>%
                mutate(
                    n_resources = as.factor(
                        sapply(
                            str_split(sources, ';'),
                            function(s){
                                s1 <- sapply(
                                    str_split(s, '_'),
                                    function(s0){tail(s0, n = 1)}
                                )
                                length(unique(s1))
                            }
                        )
                    )
                )
            
            invisible(self)
            
        }
        
    )
    
)


EnzymeSubstrateBy <- R6::R6Class(
    
    'EnzymeSubstrateBy',
    
    inherit = EnzymeSubstrateBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(by){
            
            self$by <- enquo(by)
            self$by_str <- quo_text(self$by)
            self$what_str <- `if`(
                self$by_str == 'enzyme',
                'substrate',
                'enzyme'
            )
            self$what <- parse_expr(self$what_str)
            name <- parse_expr(
                sprintf('fig_%s_by_%s', self$what_str, self$by_str)
            )
            
            super$initialize(
                name = UQ(name),
                height = 4,
                xlab_vertical = FALSE,
                theme_args = list(
                    axis.text.x = element_text(size = 14)
                ),
                width = 7
            )
            
        },
        
        
        plot = function(){
            
            DistDensHist$new(
                obj = self,
                data = self$enz_sub,
                x = counts,
                ylab = sprintf('%ss', str_to_title(self$by_str)),
                xlab = sprintf(
                    'Number of %ss\nper %s',
                    self$what_str,
                    self$by_str
                ),
                density_adjust = 2
            )
            
            invisible(self)
            
        },
            
        save = function(){
            
            super$save(print_open = FALSE)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$enz_sub <- self$enz_sub %>%
                group_by(UQ(self$by)) %>%
                mutate(counts = length(unique(UQ(self$what)))) %>%
                summarize_all(first) %>%
                ungroup()
            
            invisible(self)
            
        }
        
    )
    
)


EnzymesPerSubstrate <- R6::R6Class(
    
    'EnzymesPerSubstrate',
    
    inherit = EnzymeSubstrateBy,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            super$initialize(by = substrate, ...)
            
        }
        
    )
    
)


SubstratesPerEnzyme <- R6::R6Class(
    
    'SubstratesPerEnzyme',
    
    inherit = EnzymeSubstrateBy,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(...){
            
            super$initialize(by = enzyme, ...)
            
        }
        
    )
    
)


SubstrateNumofSites <- R6::R6Class(
    
    'SubstrateNumofSites',
    
    inherit = EnzymeSubstrateBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_substrate_numof_sites,
                height = 4,
                width = 7,
                xlab_vertical = FALSE,
                theme_args = list(
                    axis.text.x = element_text(size = 14)
                )
            )
            
        },
        
        
        plot = function(){
            
            DistDensHist$new(
                obj = self,
                data = self$enz_sub,
                x = n_sites,
                ylab = 'Substrates',
                xlab = 'Number of modification\nsites on the substrate',
                density_adjust = 2.5
            )
            
            invisible(self)
            
        },
            
        save = function(){
            
            super$save(print_open = FALSE)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$enz_sub <- self$enz_sub %>%
                group_by(
                    substrate, residue_type, residue_offset, modification
                ) %>%
                summarize_all(first) %>%
                ungroup() %>%
                group_by(substrate) %>%
                mutate(n_sites = n()) %>%
                summarize_all(first) %>%
                ungroup()
            
            invisible(self)
            
        }
        
    )
    
)


EnzymeSubstrateNumofReferences <- R6::R6Class(
    
    'EnzymeSubstrateNumofReferences',
    
    inherit = EnzymeSubstrateBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_enzyme_substrate_numof_ref,
                height = 4,
                width = 7,
                xlab_vertical = FALSE,
                theme_args = list(
                    axis.text.x = element_text(size = 14)
                ),
                keep_references = TRUE
            )
            
        },
        
        
        plot = function(){
            
            DistDensHist$new(
                obj = self,
                data = self$enz_sub,
                x = n_references,
                ylab = 'Substrate',
                xlab = paste0(
                    'Number of references\n',
                    'per enzyme-substrate interaction'
                ),
                density_adjust = 4
            )
            
            invisible(self)
            
        },
            
        save = function(){
            
            super$save(print_open = FALSE)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$enz_sub <- self$enz_sub %>%
                mutate(
                    n_references = sapply(
                        str_split(references, ';'),
                        function(s){
                            length(unique(s))
                        }
                    )
                )
            
            invisible(self)
            
        }
        
    )
    
)


EnzymeSubstratePerReference <- R6::R6Class(
    
    'EnzymeSubstratePerReference',
    
    inherit = EnzymeSubstrateBase,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            super$initialize(
                name = fig_enzyme_substrate_by_ref,
                height = 4,
                width = 7,
                xlab_vertical = FALSE,
                theme_args = list(
                    axis.text.x = element_text(size = 11)
                ),
                keep_references = TRUE
            )
            
        },
        
        
        plot = function(){
            
            DistDensHist$new(
                obj = self,
                data = self$enz_sub,
                x = n_enz_sub,
                ylab = 'References',
                xlab = paste0(
                    'Number of enzyme-substrate\n',
                    'interactions from one reference'
                ),
                density_adjust = 6
            )
            
            invisible(self)
            
        },
            
        save = function(){
            
            super$save(print_open = FALSE)
            
        }
        
    ),
    
    
    private = list(
        
        setup = function(){
            
            super$setup()
            
            self$enz_sub <- self$enz_sub %>%
                separate_rows(references, sep = ';') %>%
                group_by(references) %>%
                mutate(n_enz_sub = n()) %>%
                summarize_all(first) %>%
                ungroup()
            
            invisible(self)
            
        }
        
    )
    
)
