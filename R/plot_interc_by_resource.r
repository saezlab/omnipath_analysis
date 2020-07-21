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


IntercellNetworkByResourcePlot <- R6::R6Class(

    'IntercellNetworkByResourcePlot',

    inherit = SinglePlot,

    lock_objects = FALSE,

    public = list(

        initialize = function(
                ...
            ){

            super$initialize(
                data = self$data,
                name = fig_intercell_by_res,
                width = 6,
                height = 3,
                ...
            )

        },


        preprocess = function(...){

            invisible(self)

        },


        plot = function(...){

            self$plt <- ggplot(
                    self$data %>% filter(cnt > 0),
                    aes(x = cnt, y = resource, alpha = uni)
                ) +
                facet_grid(
                    .~obj,
                    scales = 'free_x',
                    labeller = labeller(
                        obj = c(
                            con = 'Connections',
                            receiver = 'Receivers',
                            transmitter = 'Transmitters'
                        )
                    )
                ) +
                #scale_x_continuous(expand = expansion(mult = c(.1, .1))) +
                scale_x_log10(limits=c(3, NA), expand = expansion(mult = c(.01, .07))) +
                geom_point(
                    size = 5,
                    shape = 16,
                    color = omnipath2_settings$get(palette2)[1]
                ) +
                scale_alpha_manual(
                    values = c(
                        `FALSE` = .8,
                        `TRUE` = .4
                    ),
                    labels = c(
                        `FALSE` = 'Total',
                        `TRUE` = 'Shared'
                    ),
                    guide = guide_legend(title = '')
                ) +
                ylab('Resources') +
                xlab('Number of connections and proteins')

            invisible(self)

        }

    ),


    private = list(

        setup = function(...){

            self$data <- IntercellNetworkByResource$new()$data

            invisible(self)

        }

    )

)
