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


Path <- R6::R6Class(
    
    'Path',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(
            name,
            ext = NULL,
            type = NULL,
            add_timestamp = TRUE
        ){
            
            self$name <- name
            self$type <- type
            self$ext <- `if`(is.null(ext), '', sprintf('.%s', ext))
            self$add_timestamp <- add_timestamp
            
            private$set_dir()
            private$set_fname()
            private$build_path()
            
            invisible(self)
            
        }
        
    ),
    
    private = list(
        
        get_timestamp = function(){
            
            `if`(self$add_timestamp, Sys.Date(), '')
            
        },
        
        set_dir = function(){
            
            self$dir <- case_when(
                self$type == 'figure' ~ omnipath2_settings$get(dir_figures),
                self$type == 'data' ~ omnipath2_settings$get(dir_data),
                TRUE ~ omnipath2_settings$get(dir_misc),
                NA ~ omnipath2_settings$get(dir_misc),
            )
            
            if(self$add_timestamp){
                
                self$dir <- file.path(self$dir, private$get_timestamp())
                
            }
            
            dir.create(self$dir, showWarnings = FALSE, recursive = TRUE)
            
            invisible(self)
            
        },
        
        set_fname = function(){
            
            self$fname <- sprintf(
                '%s%s%s%s',
                self$name,
                `if`(self$add_timestamp, '__', ''),
                private$get_timestamp(),
                self$ext
            )
            
            invisible(self)
            
        },
        
        build_path = function(){
            
            self$path <- file.path(self$dir, self$fname)
            
            invisible(self)
            
        }
        
    )
)


InputPath <- R6::R6Class(
    
    'InputPath',
    
    inherit = Path,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(name){
            
            super$initialize(
                name = name,
                type = 'data',
                add_timestamp = FALSE
            )
            
            invisible(self)
            
        }
        
    )
    
)


FigurePath <- R6::R6Class(
    
    'FigurePath',
    
    lock_objects = FALSE,
    
    inherit = Path,
    
    public = list(
        
        initialize = function(name){
            
            super$initialize(
                name = name,
                ext = 'pdf',
                type = 'figure'
            )
            
            invisible(self)
            
        }
        
    )
    
)
