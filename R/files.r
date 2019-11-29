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
require(jsonlite)

Files <- R6::R6Class(
    
    'Files',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(){
            
            self$files_json <- omnipath2_settings$get(files_json)
            
            self$read_files_db()
            
        },
        
        
        read_files_db = function(){
            
            self$files <- `if`(
                file.exists(self$files_json),
                jsonlite::fromJSON(self$files_json),
                list(
                    recent = list(),
                    history = list()
                )
            )
        },
        
        
        update_record = function(path){
            
            fname <- sub(
                '\\.[[:alnum:]_]+$',
                '',
                basename(path)
            )
            fname <- strsplit(fname, '__', fixed = TRUE)[[1]][1]
            
            self$read_files_db()
            
            self$files$history[[fname]] <- `if`(
                fname %in% names(self$files$history),
                unique(append(self$files$history[[fname]], path)),
                list(path)
            )
            self$files$recent[[fname]] <- path
            
            self$write_files_db()
            
        },
        
        
        write_files_db = function(){
            
            write(
                jsonlite::toJSON(self$files, pretty = TRUE),
                self$files_json
            )
            
        },
        
        
        get = function(key){
            
            self$read_files_db()
            
            return(
                `if`(
                    key %in% names(self$files$recent),
                    self$files$recent[[key]],
                    NULL
                )
            )
            
        }
        
    )
    
)
