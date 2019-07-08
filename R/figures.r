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
require(readr)
require(igraph)
require(qgraph)
require(ggraph)
require(tidygraph)
require(R6)


SinglePlot <- R6::R6Class(
    
    'SinglePlot',
    
    inherit = FigurePath,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            data,
            name,
            typeface = NA,
            width = 4,
            height = 3,
            preproc_args = list(),
            plot_args = list()
        ){
            
            super$initialize(name = name)
            
            self$data <- data
            self$width <- width
            self$height <- height
            self$preproc_args <- preproc_args
            self$plot_args <- plot_args
            
            self$main()
            
            invisible(self$plt)
        },
        
        main = function(){
            
            do.call(self$preprocess, self$preproc_args)
            do.call(self$plot, self$plot_args)
            private$set_typeface()
            private$set_theme()
            self$save()
            
        },
        
        save = function(){
            
            cairo_pdf(self$path, width = self$width, height = self$height)
            
            print(self$plt)
            
            dev.off()
            
        },
        
        preprocess = function(...){
            
            invisible(self)
            
        },
        
        plot = function(...){
            
            self$plt <- ggplot() + theme_void()
            
            invisible(self)
            
        }
        
    ),
    
    private = list(
        
        set_typeface = function(){
            
            if(is.null(self$typeface)){
                
                self$typeface <- omnipath2_settings$get(typeface)
                
            }
        },
        
        set_theme = function(){
            
            self$plt <- self$plt +
                omnipath2_settings$get(theme)() +
                theme(text = element_text(family = self$typeface))
            
        }
        
    )
    
)


PlotSeries <- R6::R6Class(
    
    'PlotSeries',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            data,
            slice_var,
            plotter,
            name,
            plot_args = list(),
            label_mapping = NULL,
            exclude_levels = NULL,
            label_levels = NULL
        ){
            
            self$data <- data
            self$name <- name
            self$slice_var <- enquo(slice_var)
            self$exclude_levels <- exclude_levels
            self$label_levels <- label_levels
            self$plotter <- plotter
            self$plot_args <- plot_args
            
            self$main()
            
            invisible(self)
            
        },
        
        
        main = function(){
            
            self$preprocess()
            self$set_levels()
            
            #self$plot()
            
        },
        
        
        preprocess = function(...){
            
            invisible(self)
            
        },
        
        
        iter_levels = function(){
            
            idx <- 0
            
            function(reset = FALSE){
                
                idx <<- `if`(reset, 1, idx + 1)
                
                `if`(
                    idx <= length(self$use_levels),
                    self$use_levels[idx],
                    NULL
                )
                
            }
            
        },
        
        
        iter_slices = function(){
            
            level_iterator <- self$iter_levels()
            
            function(reset = FALSE){
                
                level <- level_iterator(reset = reset)
                
                if(!is.null(level)){
                    
                    self$level <- level
                    self$slice_label <- `if`(
                        level %in% self$label_levels,
                        self$label_levels[[level]],
                        level
                    )
                    self$slice_name <- sprintf(
                        '%s-by-%s__%s',
                        self$name,
                        quo_name(self$slice_var),
                        level
                    )
                    self$slice <- self$data %>%
                        filter(!!self$slice_var == level)
                    
                }else{
                    
                    self$level <- NULL
                    self$slice_label <- NULL
                    self$slice_name <- NULL
                    self$slice <- NULL
                    
                }
                
            }
            
        },
        
        
        plot = function(){
            
            slice_iterator <- self$iter_slices()
            
            slice_iterator()
            
            while(!is.null(self$slice)){
                
                private$single_plot()
                
                slice_iterator()
                
            }
            
            invisible(self)
            
        },
        
        
        set_levels = function(){
            
            self$use_levels <- setdiff(
                (
                    self$data %>%
                    select(!!self$slice_var) %>%
                    c %>% unique %>% unlist
                ),
                self$exclude_levels
            )
            
        }
        
    ),
    
    private = list(
        
        single_plot = function(){
            
            do.call(
                self$plotter$new,
                c(
                    list(
                        data = self$slice,
                        name = self$slice_name
                    ),
                    self$plot_args
                )
            )
            
        }
        
    )
    
)


CategoriesPairwisePlotSeries <- R6::R6Class(
    
    'CategoriesPairwisePlotSeries',
    
    inherit = PlotSeries,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(slice_var, plotter, name, data = NULL, ...){
            
            slice_var <- enquo(slice_var)
            
            private$ensure_data(data)
            
            super$initialize(
                data = self$data,
                slice_var = !!slice_var,
                plotter = plotter,
                name = name,
                ...
            )
            
        }
        
    ),
    
    
    private = list(
        
        ensure_data = function(data){
            
            self$data <- `if`(
                is.null(data),
                IntercellCategoriesPairwise$new()$data,
                data
            )
            
            invisible(self)
            
        }
        
    )
    
)


CategorySizes <- R6::R6Class(
    
    'CategorySizes',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(data, name = NULL, ...){
            
            name <- `if`(is.null(name), 'sizes', name)
            
            super$initialize(
                data = data,
                name = name,
                ...
            )
            
        },
        
        
        preprocess = function(...){
            
            self$data <- self$data %>%
                group_by(name_cls0) %>%
                summarize_all(first) %>%
                mutate(
                    label0 = ifelse(is.na(label0), 'Total', label0)
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
                theme(
                    axis.text.x = element_text(
                        angle = 90, vjust = 0.5, hjust = 1
                    )
                ) +
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


if(FALSE){

    


    make_fr_layout <- function(g){
        # layout with qgraph
        # g is an igraph object
        el <- get.edgelist(g, names = FALSE)
        lo <- qgraph.layout.fruchtermanreingold(el, vcount = vcount(g),
                                                area = vcount(g)^2.3,
                                                repulse.rad = vcount(g)^2.1,
                                                niter = 3000)
        return(lo)
    }






    d <- suppressMessages(read_tsv('data/stats_by_resource_20190627.tsv')) %>%
        mutate(
            coverage0 = size_cls0 / omnipath0 * 100,
            coverage1 = size_cls1 / omnipath1 * 100,
            group_coverage0 = size_cls0 / size_parent0 * 100,
            group_coverage1 = size_cls1 / size_parent1 * 100
        ) %>%
        arrange(desc(size_cls0)) %>%
        mutate(
            name_cls0 = factor(name_cls0, levels = unique(name_cls0), ordered = TRUE)
        )

    pd <- d %>% filter(
        typ_cls0 == typ_cls1 &
        entity == 'protein' &
        typ_cls0 != 'misc' &
        typ_cls0 != 'small_main' &
        !is.na(name_cls0)
    )
    cd <- d %>% filter(
        typ_cls0 == typ_cls1 &
        entity == 'complex' &
        typ_cls0 != 'misc'
    )

    classes <- unique(pd$parent0)

    for(parent in classes){
        
        pd_c <- pd %>% filter(parent0 == parent) %>%
            group_by(name_cls0) %>%
            summarize_all(first) %>%
            mutate(
                label0 = ifelse(is.na(label0), 'Total', label0)
            ) %>%
            arrange(desc(size_cls0)) %>%
            mutate(
                label0 = factor(label0, levels = unique(label0), ordered = TRUE)
            )
        
        p <- ggplot(pd_c, aes(x = label0, y = size_cls0)) +
            geom_col(fill = 'black') +
            xlab('Resources') +
            ylab('Number of proteins') +
            theme_minimal() +
            theme(
                text = element_text(family = 'DINPro'),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
            ) +
            ggtitle(sprintf('Resources: %s', toupper(parent)))
        
        ggsave(
            sprintf('sizes_%s.pdf', parent),
            device = cairo_pdf,
            height = 3,
            width = 4
        )
        
    }

    edges <- d %>%
        filter(typ_cls0 == 'main' & typ_cls1 == 'main') %>%
        select(name_cls0, name_cls1, con_all, size_cls0, size_cls1) %>%
        filter(size_cls0 > 0 & size_cls1 > 0)

    vertices <- bind_rows(
            edges %>% select(name = name_cls0, size = size_cls0),
            edges %>% select(name = name_cls1, size = size_cls1)
        ) %>%
        group_by(name) %>%
        summarize_all(first)

    edges <- edges %>% select(name_cls0, name_cls1, con_all)

    g <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE)
    V(g)$size <- V(g)$size / max(V(g)$size) * 50 + 30
    E(g)$width <- E(g)$con_all / max(E(g)$con_all) * 20 + 3

    g <- simplify(g, edge.attr.comb = first, remove.loops = FALSE)

    lo <- make_fr_layout(g)

    cairo_pdf('figures/intercell_graph.pdf', height = 5, width = 5)

        plot(
            g,
            layout = lo,
            vertex.frame.color = NA,
            vertex.color = '#FDC70F',
            edge.color = '#4C4B6B33',
            vertex.label.family = 'DINPro',
            vertex.label.color = '#4C4B6B',
            edge.width = E(g)$width
        )

    dev.off()

}
