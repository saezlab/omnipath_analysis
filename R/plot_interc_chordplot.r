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
require(circlize)


IntercellChordplot <- R6::R6Class(

    'IntercellChordplot',

    inherit = SinglePlot,

    lock_objects = FALSE,

    public = list(

        initialize = function(input_param = NULL, ...){

            self$input_param <- `if`(
                is.null(input_param),
                c('omnipath', 'undirected', 'any_effect'),
                input_param
            )

            super$initialize(
                data = self$data,
                name = fig_interc_chordplot,
                fname_param = self$input_param,
                height = 4,
                width = 4
            )

            invisible(self)

        },


        main = function(){

            private$setup()
            self$plot()

            invisible(self)

        },


        plot = function(...){

            grid.col <- omnipath2_settings$get(intercell_colors)
            grid.col <- grid.col[names(grid.col) %in% self$data$from] %>%
                unlist()

            self$open_device()
            self$plt <- chordDiagram(
                self$data,
                grid.col = grid.col,
                transparency = .3,
                self.link = 2
            )
            dev.off()
            circos.clear()
            private$ready()

            invisible(self)

        },


        post_plot = function(){

            invisible(self)

        }

    ),


    private = list(

        setup = function(){

            self$data <- IntercellNetworkSummary$new(
                    input_param = self$input_param
                )$data %>%
                rename(
                    from = category_a,
                    to = category_b,
                    value = connections
                )

            invisible(self)

        }

    )

)