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
            
            invisible(self)
        },
        
        main = function(){
            
            private$setup()
            do.call(self$preprocess, self$preproc_args)
            do.call(self$plot, self$plot_args)
            private$post_plot()
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
        
        setup = function(){
            
            invisible(self)
            
        },
        
        
        post_plot = function(){
            
            private$set_typeface()
            private$set_theme()
            
            invisible(self)
            
        },
        
        
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
