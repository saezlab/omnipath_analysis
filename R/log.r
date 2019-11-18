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


Logger <- R6::R6Class(
    
    'Logger',
    
    lock_objects = FALSE,
    
    public = list(
        
        initialize = function(level = NULL){
            
            self$logdir <- omnipath2_settings$get('dir_log')
            self$session_id <- omnipath2_settings$get('session_id')
            self$level <- `if`(
                is.null(level),
                omnipath2_settings$get('loglevel'),
                level
            )
            
            private$init_log()
            
            invisible(self)
            
        },
        
        
        msg = function(msg, label = NULL, level = 0){
            
            if(level <= self$level){
                
                message <- sprintf(
                    '[%s] [%s] %s',
                    private$get_timestamp(),
                    `if`(is.null(label), 'omnipath2', label),
                    msg
                )
            }
            
            write(message, self$conn)
            
            invisible(self)
            
        },
        
        
        finalize = function(){
            
            private$close()
            
        }
        
    ),
        
    
    private = list(
        
        init_log = function(){
            
            private$set_path()
            private$create_dir()
            private$open()
            private$welcome()
            
            invisible(self)
            
        },
        
        
        set_path = function(){
            
            
            self$logfile <- file.path(
                self$logdir,
                sprintf('omnipath2-%s.log', self$session_id)
            )
            
            invisible(self)
            
        },
        
        
        create_dir = function(){
            
            dir.create(self$logdir, showWarnings = FALSE, recursive = TRUE)
            
            invisible(self)
            
        },
        
        
        open = function(){
            
            self$conn <- file(self$logfile, open = 'at')
            
            invisible(self)
            
        },
        
        
        welcome = function(){
            
            self$msg('OmniPath2 plotting package reloaded.')
            
            invisible(self)
            
        },
        
        
        get_timestamp = function(){
            
            format(Sys.time(), '%Y-%m-%d %H:%M:%S')
            
        },
        
        
        close = function(){
            
            close(self$conn)
            
            invisible(self)
            
        }
        
    )
    
)
