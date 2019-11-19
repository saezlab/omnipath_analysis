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
            fname_param = list(),
            ext = NULL,
            type = NULL,
            add_timestamp = NA
        ){
            
            self$name <- name
            self$type <- type
            self$ext <- `if`(is.null(ext), '', ext)
            self$add_timestamp <- add_timestamp
            self$fname_param <- fname_param
            
            private$set_dir()
            private$set_fname()
            private$build_path()
            
            invisible(self)
            
        }
        
    ),
    
    private = list(
        
        get_timestamp = function(){
            
            format(Sys.Date(), omnipath2_settings$get(timestamp_format))
            
        },
        
        set_dir = function(){
            
            self$dir <- case_when(
                self$type == 'figure' ~
                    omnipath2_settings$get(dir_figures),
                self$type %in% c('table', 'data') ~
                    omnipath2_settings$get(dir_tables),
                TRUE ~
                    omnipath2_settings$get(dir_misc),
                NA ~
                    omnipath2_settings$get(dir_misc),
            )
            
            self$dir <- `if`(
                `||`(
                    self$add_timestamp,
                    omnipath2_settings$get(timestamp_dirs)
                ),
                file.path(self$dir, private$get_timestamp()),
                self$dir
            )
            
            dir.create(self$dir, showWarnings = FALSE, recursive = TRUE)
            
            invisible(self)
            
        },
        
        timestamp_fname = function(fname){
            
            fname <- `if`(
                is_quosure(fname),
                quo_get_expr(fname),
                fname
            )
            
            fname <- `if`(
                `||`(
                    self$add_timestamp,
                    omnipath2_settings$get(timestamp_files)
                ),
                sprintf('%s__%s', fname, private$get_timestamp()),
                fname
            )
            
            return(fname)
            
        },
        
        
        set_fname_key = function(){
            
            self$fname_key <- omnipath2_settings$get(UQ(self$name))
            
            invisible(self)
            
        },
        
        
        insert_param = function(){
            
            self$fname_key <- do.call(
                sprintf,
                c(
                    list(self$fname_key),
                    self$fname_param
                )
            )
            
            invisible(self)
            
        },
        
        
        set_fname = function(){
            
            private$set_fname_key()
            private$insert_param()
            
            self$fname <- sprintf(
                '%s.%s',
                private$timestamp_fname(self$fname_key),
                self$ext
            )
            
            invisible(self)
            
        },
        
        build_path = function(){
            
            self$path <- file.path(self$dir, self$fname)
            
            invisible(self)
            
        },
        
        
        ready = function(){
            
            omnipath2_files$update_record(self$path)
            
            invisible(self)
            
        }
        
    )
)


InputPath <- R6::R6Class(
    
    'InputPath',
    
    inherit = Path,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(name, ...){
            
            super$initialize(
                name = enquo(name),
                fname_param = list(...),
                type = 'data',
                add_timestamp = FALSE,
                ext = 'tsv'
            )
            
            invisible(self)
            
        }
        
    ),
    
    
    private = list(
        
        set_dir = function(){
            
            return(self)
            
        },
        
        
        set_fname = function(){
            
            private$set_fname_key()
            private$insert_param()
            
            self$path <- omnipath2_files$get(self$fname_key)
            self$dir <- dirname(self$path)
            self$fname <- basename(self$path)
            
            invisible(self)
            
        },
        
        
        build_path = function(){
            
            invisible(self)
            
        }
        
    )
    
)


TablePath <- R6::R6Class(
    
    'TablePath',
    
    inherit = Path,
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(name, fname_param = list()){
            
            super$initialize(
                name = name,
                fname_param = fname_param,
                type = 'table',
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
        
        initialize = function(name, fname_param = list()){
            
            super$initialize(
                name = name,
                fname_param = fname_param,
                ext = 'pdf',
                type = 'figure'
            )
            
            op2log(
                sprintf('Figure path set: `%s`.', self$path),
                label = class(self)[1]
            )
            
            invisible(self)
            
        }
        
    )
    
)
