proc ::DBTools::setGeometry\
  {selections molids atoms geometry} {

      ## This is a bit clunky, but it is the minimal
      ## number of transformations required to get the
      ## geometry about the linkage correct

      if {[llength $geometry] != 6} {
        dbtCon -error "Geometry list format is: chi1 r1 theta1 chi2 theta2 chi3"
        return -code error
      }

      if {[llength $atoms] != 6} {
        dbtCon -error "six atoms are required for link definition"
        return -code error
      }

      set nmolids [llength $molids] 

      if {$nmolids == 2} {
 
        lassign $atoms a1 a2 a3 a4 a5 a6
        lassign $molids m1 m2
        set midsatoms [list $a1 $m1 $a2 $m1 $a3 $m1\
          $a4 $m2 $a5 $m2 $a6 $m2]
      
      } elseif {$nmolids == 6} {
    
        set midsatoms {}
        foreach a $atoms m $molids { lappend midsatoms $a $m} 
      
      } else {
          dbtCon -error "two or six molids are required for geometry specification"
          return -code error
      }

      ## Get selections and geometries 
      set sels [getSel selections $midsatoms] 
      lassign $sels sel1 sel2 sel3 sel4 sel5 sel6 
      lassign $geometry chi1 r theta1 chi2 theta2 chi3

      ## Set the geometries 
      if {$r != ""}      {setBond  [list $sel3 $sel4] $r}
      if {$chi1 != ""}   {setDihed [list $sel1 $sel2 $sel3 $sel4] $chi1}
      if {$theta1 != ""} {setAngle [list $sel2 $sel3 $sel4] $theta1}
      if {$chi2 != ""}   {setDihed [list $sel2 $sel3 $sel4 $sel5] $chi2}
      if {$theta2 != ""} {setAngle [list $sel3 $sel4 $sel5] $theta2}
      if {$chi3 != ""}   {setDihed [list $sel3 $sel4 $sel5 $sel6] $chi3}
}


proc ::DBTools::getSel {selections midsatoms} {

  ## Return set of selections
  ## used to set the geometry
  ## between the specified fragments

  set retval {}
  foreach {atom molid} $midsatoms {
    if {$atom == "" || $molid == ""} {
      lappend retval ""
    } else { 
      lappend retval [__getSel selections $molid $atom]
    }
  }

  return $retval
}

proc ::DBTools::__getSel {selections mol atom} {

  upvar #1 $selections s

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

proc ::DBTools::setBond {sels r} {
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

proc ::DBTools::setAngle {sels theta} {
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

proc ::DBTools::setDihed {sels chi} {
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
