#!/usr/bin/tclsh

# DBTools. Tools for building dendrimers 

## None of this will work without topotools.
package require topotools

## This package makes extensive use of some Tcl 8.5 features.
package require Tcl 8.5

namespace eval ::DBTools:: {

    namespace export DBTools

    variable version 0.1

    ## Store interface variables
    variable sys
    array unset sys

    ## Show the title and citation info?
    variable showtitle 0

    ## Locations of additional data files
    variable datadir $env(DBTOOLSDIR)

    ## Return error codes
    set sys(OK) 0
    set sys(ERROR) -1
}

# +----------------+
# | Global command |
# +----------------+

proc DBTools { args } {
    eval ::DBTools::DBTools $args
}

# +---------+
# | Startup |
# +---------+

interp alias {} dbt {} ::DBTools::DBTools
package provide dbt $::DBTools::version

if {$::DBTools::showtitle} {::DBTools::title}

## Load other cgtools files
source [file join $env(DBTOOLSDIR) dbtools_geometry.tcl]
