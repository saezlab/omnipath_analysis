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
require(igraph)
require(qgraph)
require(ggraph)
require(tidygraph)
require(R6)


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
