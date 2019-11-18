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

require(R6)


omnipath2_workflow <<- list(
    
    # These "tests" show if we can provide the arguments
    # to the workflow elements using Task and Param
    # classes.
    # Here we just read 2 types of data tables,
    # a number of variants each.
    test1 = Task$new(
        method = ResourceByEntity,
        name = NULL,
        Param$new('omnipath', 'curated')
    ),
    
    test2 = Task$new(
        method = IntercellCategoriesPairwise,
        name = NULL,
        ProductParam$new(
            list('curated', 'omnipath'),
            list('proteins'),
            list(
                'all-class-levels',
                'above_main-main-misc-small_main'
            )
        )
    ),
    
    test3 = Task$new(
        method = IntercellCategoriesPairwise,
        name = NULL,
        Param$new(NULL),
        Param$new(
            'curated', 'omnipath'
        ),
        Param$new('proteins'),
        Param$new(
            'all-class-levels',
            'above_main-main-misc-small_main'
        )
    )
    
)
