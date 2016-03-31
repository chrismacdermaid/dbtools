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
    array unset sys *

    ## Show the title and citation info?
    variable showtitle 0

    ## Locations of additional data files
    variable datadir $env(DBTOOLSDIR)

    ## Return error codes
    set sys(OK) 0
    set sys(ERROR) -1

    ## Properties to copy from templates
    variable cpylist {name type mass charge radius element \
      resname resid chain segname x y z}

    ## Array containing the binary tree connectivities
    variable tree
    array unset tree *

    ## Dendrimer Topologies 
    variable topo
    set topo {}
}

## A wrapper around cgCon
proc DBTools::dbtCon {flag str} {
    switch -- $flag {
        "-error" -
        "-err"  { vmdcon -err  "DBTOOLS> $str" }
        "-warn" { vmdcon -warn "DBTOOLS> $str" }
        "-info" { vmdcon -info "DBTOOLS> $str" }
        default { vmdcon "DBTOOLS> $str"}
    }
}

# +----------------+
# | Global command |
# +----------------+

proc DBTools { args } {
    eval ::DBTools::DBTools $args
}

proc DBTools::DBTools args {return -code ok}

# +---------+
# | Startup |
# +---------+

interp alias {} dbt {} ::DBTools::DBTools
package provide dbt $::DBTools::version

if {$::DBTools::showtitle} {::DBTools::title}

## Load other cgtools files
#set ::env(DBTOOLSDIR) [pwd] ;# for testing
source [file join $env(DBTOOLSDIR) dbtools_build.tcl]
source [file join $env(DBTOOLSDIR) dbtools_geometry.tcl]
source [file join $env(DBTOOLSDIR) dbtools_tree.tcl]
