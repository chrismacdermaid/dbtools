
proc ::DBTools::buildDendrimer {dname gen\
  {coreid ""} {repeatid ""} {termid ""}} {

  ## Build a nth generation dendrimer
  ## having name "dname" defined in the 
  ## dendrimer topology

  variable topo 
  variable tree

  ## Load Fragments if necessary 
  set names [dict get $topo $dname name] 
  if {$coreid == ""} {
   set prefix [dict get $names core]
   set coreid [safeNewMol $prefix.pdb]
  }

  if {$repeatid == ""} {
   set prefix [dict get $names repeat]
   set repeatid [safeNewMol $prefix.pdb]
  }

  if {$termid == ""} {
   set prefix [dict get $names term]
   set termid [safeNewMol $prefix.pdb]
  }

  ## Create a dummy molecule to store temprary files 
  set dumid [safeNewMol [molinfo $repeatid get filename]]

  ## Create some selections for each of the fragments
  set sel_core [atomselect $coreid "all"]
  set sel_repeat [atomselect $repeatid "all"]
  set sel_dum [atomselect $dumid "all"] 
  if {$termid != ""} {
    set sel_term [atomselect $termid "all"]
  }

  ## Coordinates, topologies 
  set coords {}; set top {}; set bondlist {} 

  ## Root id of the binary tree 
  set rootid [lindex $tree($gen) 0 0] 

  ## Center the core at the origin
  $sel_core moveby [vecinvert [measure center $sel_core]]
  set coord_arr($rootid) [$sel_core get {x y z}]
  
  ## update the data
  lappend top [$sel_core get $cpylist]
  lappend bondlist [topo -sel $sel_core getbondlist] 

  ## Attach the core to the repeats 
  set links [dict get $topo $dname link core-repeat]
  set geometry [dict get $topo $dname geometry core-repeat]
  
  ## Get the nodes with depth = 1 
  set paths  
  foreach p $paths l $links g $geometry {
   set id0 [lindex $p 0]; ## Parent 
   set id1 [lindex $p 1]; ## Left leaf
  
   ## Load the coordinates of the antecendent
   ## into the dummy mol
   $sel_dum set {x y z} $coord_arr($id0)
  
   ## Set the Geometry about the first linkage 
   setGeometry selections $l $g

   ## Update coordiates, topologies and properties
   set coord_arr($id1) [$sel_repeat get {x y z}]
   $sel_repeat set resid "$root$id1"
   lappend top [$sel_repeat get $cpylist] 
   lappend bondlist [topo -sel $sel_repeat getbondlist] 
  }
  
  ## Generate the coordinate list
  foreach {key value} [array get coord_arr *] {
    lappend coords $value
  }

  ## Collapse the coordinate and property arrays
  set coords [join $coords]
  set top [join $top]

  ## Merge the fragments 
  set newmol [mergeFragments $coords $top 

  ## Fix bonds

  return $newmol

}

proc ::DBTools::mergeFragments {coords top} {

    set coords [join $coords]
    set top [join $top]

    ## Create a new mol and dump the coordinates
    ## and atom properties into it
    set newmol [mol new atoms [llength $coords]]
    animate dup $newmol
    set sel_new [atomselect $newmol "all"]
    $sel_new set {x y z} $coords
    $sel_new set $cpylist $top

    ## Update the bond list and create bonds between
    ## fragments 
    set newbl {}
    set offset [llength [lsort\
        -unique -integer -increasing [join [lindex $bondlist 0]]]]
    lappend newbl {*}[lindex $bondlist 0]

    foreach l [lrange $bondlist 1 end] {
      ## Collapse the list and get the unique
      ## indice, count them
        foreach b $l {
          lassign $b bid1 bid2
            incr bid1 $offset
            incr bid2 $offset
            lappend newbl [list $bid1 $bid2]
        }

      incr offset [llength [lsort\
        -unique -integer -increasing [join $l]]]
    }

    ## Apply the bond list
    topo -molid $newmol setbondlist $newbl

    ## Cleanup
    $sel_new delete

    return $newmol
}

proc ::DBTools::safeMolNew {fname} {

  if { [catch { mol new $fname } retval } {
    dbtCon -error $retval; return -code error
  }

  return $retval 
}
