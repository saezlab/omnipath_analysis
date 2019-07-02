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
            plt,
            name,
            typeface = NA,
            width = 4,
            height = 3
        ){
            
            super$initialize(name = name)
            
            self$plt <- plt
            self$width <- width
            self$height <- height
            
            self$main()
            
            invisible(self$plt)
        },
        
        main = function(
            
        ){
            
            self$set_typeface()
            self$set_theme()
            self$save()
            
        },
        
        save = function(){
            
            cairo_pdf(self$path, width = self$width, height = self$height)
            
            print(self$plt)
            
            dev.off()
            
        }
        
    ),
    
    private = list(
        
        set_typeface = function(){
            
            if(is.na(self$typeface)){
                
                self$typeface <- settings$get(typeface)
                
            }
        },
        
        set_theme = function(){
            
            self$plt <- self$plt +
                omnipath2_settings$get(theme)() +
                theme(text = element_text(family = self$typeface))
            
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
