## Build a nth generation of dendrimer

proc buildDendrimer {gen coreid repeatid {termid none}} {

    catch {array unset selections}
    array set selections {}

    ## Create a dummy molecule
    set dumid [mol new [molinfo $repeatid get filename]] 
    
    ## Make some selections for the pieces
    set sel_core [atomselect $coreid "all"]
    set sel_repeat [atomselect $repeatid "all"]
    set sel_dum [atomselect $dumid "all"] 
    if {$termid != "none"} {
      set sel_term [atomselect $termid "all"]
    }

    ## Coordinates, topologies 
    set coords {}
    set top {}
    set bondlist {} 

    set cpylist {name type mass charge radius element \
      resname resid chain segname}

    ## Center the core at the origin
    $sel_core moveby [vecinvert [measure center $sel_core]]
    lappend coords [$sel_core get {x y z}]
    lappend top [$sel_core get $cpylist]
    lappend bondlist [topo -sel $sel_core getbondlist] 

    set core1_repeat [list [list C2 $coreid] [list C1 $coreid]\
          [list N1 $coreid] [list C1 $repeatid]\
          [list C2 $repeatid] [list C5 $repeatid]]
    setGeometry selections $core1_repeat {175.0 1.35 125.7 -170 115.3 -76.0} 
    lappend coords [$sel_repeat get {x y z}]
    lappend top [$sel_repeat get $cpylist]
    lappend bondlist [topo -sel $sel_repeat getbondlist] 

    set core2_repeat [list [list C1 $coreid] [list C2 $coreid]\
          [list N2 $coreid] [list C1 $repeatid]\
          [list C2 $repeatid] [list C5 $repeatid]]
    setGeometry selections $core2_repeat {175.0 1.35 125.7 -170 115.3 -76.0} 
    lappend coords [$sel_repeat get {x y z}]
    lappend top [$sel_repeat get $cpylist]
    lappend bondlist [topo -sel $sel_repeat getbondlist] 

    ## Load in the binary trees that describe the dendrimer
    ## connectivity
    source trees.dat

    ## Now the tricky part. We need to traverse a binary tree
    ## and use the correct geometries 
    catch {array unset coord_arr}
    array set coord_arr {}
 
    ## Specified generation should be one lower for these
    ## since were building each hemisphere separately
    incr gen -1

    if {1} { 

        set repeat1_repeat [list [list C2 $dumid] [list C3 $dumid]\
                           [list N1 $dumid] [list C1 $repeatid]\
                           [list C2 $repeatid] [list C5 $repeatid]]

        set repeat2_repeat [list [list C2 $dumid] [list C4 $dumid]\
                           [list N2 $dumid] [list C1 $repeatid]\
                           [list C2 $repeatid] [list C5 $repeatid]]

      foreach root {1 2} {  

        ## Get the root of the tree's ID 
        set coord_arr([lindex $arr($gen) 0 0]) [lindex $coords $root]

        for {set i 1} {$i <= $gen} {incr i} {
          set paths [lsort -unique -integer -index $i $arr($gen)]
          foreach {a b} $paths {
            set id0 [lindex $a $i-1]; ## Parent 
            set id1 [lindex $a $i]; ## Left leaf
            set id2 [lindex $b $i]; ## Right leaf 

            ## If the coordinates don't exist
            ## we need to generate them using the antecendent
            if {![info exists coord_arr($id1)]} {
              ## Load the coordinates of the antecendent
              ## into the dummy mol
              $sel_dum set {x y z} $coord_arr($id0)
              ## Set the link-geometry
              setGeometry selections $repeat1_repeat {175.0 1.35 125.7 -170 115.3 -76.0} 
              set coord_arr($id1) [$sel_repeat get {x y z}]
              $sel_repeat set resid "$root$id1"
              lappend top [$sel_repeat get $cpylist] 
              lappend bondlist [topo -sel $sel_repeat getbondlist] 
            }    

            if {![info exists coord_arr($id2)]} {
              ## Load the coordinates of the antecendent
              ## into the dummy mol
              $sel_dum set {x y z} $coord_arr($id0)
              ## Set the link-geometry
              setGeometry selections $repeat2_repeat {175.0 1.35 125.7 -170 115.3 -76.0} 
              set coord_arr($id2) [$sel_repeat get {x y z}]
              $sel_repeat set resid "$root$id2"
              lappend top [$sel_repeat get $cpylist] 
              lappend bondlist [topo -sel $sel_repeat getbondlist] 
            }    
          }
        } 
    
        ## Generate the coordinate list
        foreach {key value} [array get coord_arr *] {
          ## Don't double add the roots
          if {$key == [lindex $arr($gen) 0 0]} {continue}
          lappend coords $value
        }

         ## Clear the coordinate array 
         array unset coord_arr * 
    }
  }
    mol delete $dumid

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

    topo -molid $newmol setbondlist $newbl



    ## Cleanup
    foreach {key value} [array get selections] {
        catch {$value delete}
        unset selections($key)
    }
}

## Helpers

proc setGeometry {selections midsatoms geometry} {

      ## This is a bit clunky, but it is the minimal
      ## number of transformations required to get the
      ## geometry about the linkage correct

      if {[llength $geometry] != 6} {
        puts "Geometry list format is: chi1 r1 theta1 chi2 theta2 chi3"
        return -code error
      }

      set sels [getSel selections $midsatoms] 
      lassign $sels sel1 sel2 sel3 sel4 sel5 sel6 
      
      lassign $geometry chi1 r theta1 chi2 theta2 chi3

      if {$r != ""}      {setBond  [list $sel3 $sel4] $r}
      if {$chi1 != ""}   {setDihed [list $sel1 $sel2 $sel3 $sel4] $chi1}
      if {$theta1 != ""} {setAngle [list $sel2 $sel3 $sel4] $theta1}
      if {$chi2 != ""}   {setDihed [list $sel2 $sel3 $sel4 $sel5] $chi2}
      if {$theta2 != ""} {setAngle [list $sel3 $sel4 $sel5] $theta2}
      if {$chi3 != ""}   {setDihed [list $sel3 $sel4 $sel5 $sel6] $chi3}
}


proc getSel {selections midsatoms} {

  ## Return set of selections
  ## used to set the geometry
  ## between the specified fragments

  set retval {}
  foreach x $midsatoms {
    lassign $x atom molid
    lappend retval [__getSel selections $molid $atom] 
  }

  return $retval
}

proc __getSel {selections mol atom} {
  upvar #1  $selections s

  ## Check if a selection exists and 
  ## return handle if it does; otherwise
  ## create it and return handle

  if {[info exists s([list $mol $atom])]} {
    set retval $s([list $mol $atom])
  } else {
    set retval [atomselect $mol "name $atom"]; $retval global 
    set s([list $mol $atom]) $retval
  }

  return $retval
}

proc setBond {sels r} {
  ## Set the bond length between two atoms

  lassign $sels sel1 sel2

  ## Get the coordinates of the atoms
  set xyz1 [join [$sel1 get {x y z}]]
  set xyz2 [join [$sel2 get {x y z}]]

  set v [vecsub $xyz1 $xyz2]

  set frag2 [atomselect [$sel2 molid]\
    "same fragment as index [$sel2 get index]"]

  $frag2 moveby $v
  $frag2 moveby [vecscale\
    [vecinvert [vecnorm $v]] $r] 
  $frag2 delete

  return -code ok 
}

proc setAngle {sels theta} {
  ## Set the angle for three atoms

  lassign $sels sel1 sel2 sel3

  ## Atom indices
  set idx1 [$sel1 get index] 
  set idx2 [$sel2 get index] 
  set idx3 [$sel3 get index] 

  ## Get the coordinates of the atoms
  set xyz1 [join [$sel1 get {x y z}]]
  set xyz2 [join [$sel2 get {x y z}]]
  set xyz3 [join [$sel3 get {x y z}]]

  ## Measure the current angle value
  set theta0 [measure angle [list\
  [list $idx1 [$sel1 molid]]\
  [list $idx2 [$sel2 molid]]\
  [list $idx3 [$sel3 molid]]]]

  ## Calculate the offset
  set dtheta [expr {$theta0 - $theta}]

  ## Rotate to set the desired angle  
  set r [trans angle $xyz1 $xyz2 $xyz3 $dtheta deg]

  set frag3 [atomselect [$sel3 molid]\
    "same fragment as index [$sel3 get index]"] 
  $frag3 move $r
  $frag3 delete

  return -code ok 
}

proc setDihed {sels chi} {
  ## Set the dihedral for four atoms

  lassign $sels sel1 sel2 sel3 sel4

  ## Atom indices
  set idx1 [$sel1 get index] 
  set idx2 [$sel2 get index] 
  set idx3 [$sel3 get index] 
  set idx4 [$sel4 get index] 

  ## Get the coordinates of the atoms
  set xyz1 [join [$sel1 get {x y z}]]
  set xyz2 [join [$sel2 get {x y z}]]
  set xyz3 [join [$sel3 get {x y z}]]
  set xyz4 [join [$sel4 get {x y z}]]

  ## Measure the current dihedral value
  set chi0 [measure dihed [list\
  [list $idx1 [$sel1 molid]]\
  [list $idx2 [$sel2 molid]]\
  [list $idx3 [$sel3 molid]]\
  [list $idx4 [$sel4 molid]]]]

  ## Calculate the offset
  set dchi [expr {$chi - $chi0}]

  ## Rotate about atom2-atom3
  set r [trans bond $xyz2 $xyz3 $dchi deg]

  set frag4 [atomselect [$sel4 molid]\
   "same fragment as index [$sel4 get index]"]
  $frag4 move $r
  $frag4 delete

  return -code ok

}
