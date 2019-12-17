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
require(gridExtra)
require(R6)


DistDensHist <- R6::R6Class(
    
    'DistDensHist',
    
    lock_objects = FALSE,
    
    public = list(

        initialize = function(
            obj,
            data,
            x,
            ylab,
            xlab
        ){
            
            self$obj <- obj
            self$data <- data
            self$x <- enquo(x)
            self$ylab <- ylab
            self$xlab <- xlab
            
            self$main()
            
            invisible(self)
            
        },
        
        main = function(){
            
            self$obj$open_device()
            plot.new()
            
            self$obj$plt <- ggplot(self$data) +
                stat_bin(
                        aes(x = !!self$x, y = cumsum(..count..)),
                        geom = 'step',
                        binwidth = .01
                ) +
                xlab(self$xlab) +
                ylab(self$ylab) +
                scale_x_log10(expand = c(.01, 0)) +
                annotation_logticks(sides = 'b')
            
            self$obj$post_plot()
            self$obj$dist <- self$obj$plt
            
            self$obj$plt <- ggplot(self$data) +
                geom_density(aes(x = !!self$x, y = ..count..), adjust = 2) +
                ylab(self$ylab) +
                xlab(self$xlab) +
                scale_x_log10(expand = c(.01, 0)) +
                annotation_logticks(sides = 'b')
            
            self$obj$post_plot()
            self$obj$dens <- self$obj$plt
            
            self$obj$plt <- ggplot(self$data) +
                geom_histogram(
                    aes(x = !!self$x, y = ..count..),
                    fill = 'black',
                    binwidth = .25
                ) +
                ylab(self$ylab) +
                xlab(self$xlab) +
                scale_x_log10(expand = c(.01, 0)) +
                annotation_logticks(sides = 'b')
            
            self$obj$post_plot()
            self$obj$hist <- self$obj$plt
            
            self$obj$plt <- grid.arrange(
                self$obj$dist,
                self$obj$dens,
                self$obj$hist,
                nrow = 1
            )
            
            invisible(self)
        
        }
        
    )
    
)
