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

.onLoad <- function(libname, pkgname){
    
    if(!exists('omnipath2_session')){
        
        omnipath2_session <<- gen_session_id()
        
    }
    
    omnipath2_settings <<- Settings$new()
    omnipath2_settings$set(session_id = omnipath2_session)
    do.call(options, omnipath2_settings$get(console_settings))
    omnipath2_log <<- Logger$new()
    
    omnipath2_files <<- Files$new()
    
    op2log <<- function(...){
        
        omnipath2_log$msg(...)
        
    }
    
}
