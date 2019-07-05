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


Settings <- R6::R6Class(
    
    'Settings',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            ...
        ){
            
            self$reset(...)
            
            invisible(self)
            
        },
        
        get = function(key){
            
            key <- quo_name(enquo(key))
            
            self$settings[[key]]
            
        },
        
        set = function(...){
            
            self$settings <- modifyList(self$settings, list(...))
            
        },
        
        reset = function(...){
            
            self$settings <- modifyList(get_default_param(), list(...))
            
        },
        
        print = function(){
            
            cat('Current settings:\n')
            
            for(key in names(self$settings)){
                
                cat(sprintf('\t%s\t\t%s\n', key, self$settings[[key]]))
                
            }
            
        }
        
    )
    
)
