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
require(igraph)
require(qgraph)
require(ggraph)
require(tidygraph)
require(R6)


GraphPlot <- R6::R6Class(
    
    'GraphPlot',
    
    inherit = SinglePlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(data, name, layout = 'fr', ...){
            
            assign('layout_types', c('fr', 'circle'), envir = private$static)
            
            self$layout <- layout
            
            super$initialize(
                data = data,
                name = name,
                ...
            )
            
            invisible(self)
            
        },
        
        
        set_fr_layout = function(){
            
            self$layout <- self$fr_layout(self$graph)
            
            invisible(self)
            
        },
        
        
        set_circle_layout = function(){
            
            self$layout <- layout_in_circle(self$graph)
            
            invisible(self)
            
        },
        
        
        #' Layout with `qgraph`.
        #' `g` is an `igraph` object.
        fr_layout = function(g){
            
            el <- get.edgelist(g, names = FALSE)
            lo <- qgraph.layout.fruchtermanreingold(
                el,
                vcount = vcount(g),
                area = vcount(g)^2.3,
                repulse.rad = vcount(g)^2.1,
                niter = 3000
            )
            
            return(lo)
            
        },
        
        
        preprocess = function(...){
            
            invisible(self)
            
        },
        
        
        plot = function(...){
            
            self$plot_args['layout'] <- self$layout
            
            cairo_pdf(self$path, width = self$width, height = self$height)
            
            do.call(
                plot,
                c(list(self$graph), self$plot_args)
            )
            
            dev.off()
            
            invisible(self)
            
        },
        
        
        save = function(){
            
            invisible(self)
            
        },
        
    ),
    
    
    private = list(
        
        static = new.env(),
        
        
        setup = function(){
            
            private$ensure_graph()
            self$plot_args <- modifyList(
                omnipath2_settings$get(graph_plot_defaults),
                self$plot_args
            )
            private$set_layout()
            private$set_typeface()
            self$plot_args <- modifylist(
                list(
                    vertex.label.family = self$typeface,
                    edge.label.family = self$typeface,
                    layout = self$layout
                ),
                self$plot_args
            )
            
            invisible(self)
            
        },
        
        
        ensure_graph = function(){
            
            self$graph <- `if`(
                class(self$data) == 'igraph',
                self$data,
                private$build_graph()
            )
            
            invisible(self)
            
        },
        
        
        build_graph = function(){
            
            make_empty_graph()
            
        },
        
        
        post_plot = function(){
            
            invisible(self)
            
        },
        
        
        set_layout = function(){
            
            self$layout <- `if`(
                (
                    !is.null(dim(self$layout)) &
                    dim(self$layout)[1] == vcount(self$graph)
                ),
                self$layout,
                self[[
                    sprintf(
                        'set_%s_layout',
                        `if`(
                            (
                                is.character(self$layout) &
                                self$layout %in% get(
                                    'layout_types',
                                    envir = private$static
                                )
                            ),
                            self$layout,
                            omnipath2_settings$get(graph_layout_default)
                        )
                    )
                ]]()
            )
            
            invisible(self)
            
        }
        
    )
    
)


ConnectionGraph <- R6::R6Class(
    
    'ConnectionGraph',
    
    inherit = GraphPlot,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(data = NULL, layout = 'fr', ...){
            
            private$ensure_data(data)
            
            super$initialize(
                data = data,
                name = 'intercell-classes-graph',
                layout = layout,
                ...
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        ensure_data = function(data = NULL){
            
            self$data <- `if`(
                is.null(data),
                IntercellCategoriesPairwise$new()$data,
                data
            )
            
            invisible(self)
            
        },
        
        
        build_graph = function(){
            
            edges <- self$data %>%
                filter(typ_cls0 == 'main' & typ_cls1 == 'main') %>%
                select(
                    name_cls0, name_cls1, con_all, size_cls0, size_cls1
                ) %>%
                filter(size_cls0 > 0 & size_cls1 > 0)
            
            vertices <- bind_rows(
                    edges %>% select(name = name_cls0, size = size_cls0),
                    edges %>% select(name = name_cls1, size = size_cls1)
                ) %>%
                group_by(name) %>%
                summarize_all(first)
            
            edges <- edges %>% select(name_cls0, name_cls1, con_all)
            
            g <- graph_from_data_frame(
                edges,
                vertices = vertices,
                directed = FALSE
            )
            
            V(g)$size <- V(g)$size / max(V(g)$size) * 50 + 30
            E(g)$width <- E(g)$con_all / max(E(g)$con_all) * 20 + 3
            
            g <- simplify(g, edge.attr.comb = first, remove.loops = FALSE)
            
            return(g)
            
        }
        
    )
    
)
