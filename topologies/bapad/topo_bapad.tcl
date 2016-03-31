if {![namespace exists ::DBTools]} {
  namespace eval ::DBTools:: {}
  set testflag 1
}

## Define core, repeat and terminal resnames 
## and their linkages

#proc ::DBTools::topo_bapad {} {
#
#  variable topo
#
#  ## Fragment Names 
#  set topo([list core bapad])   EDH 
#  set topo([list repeat bapad]) AMP 
#  set topo([list term bapad])   AHH 
#
#  ## Names of atoms invloved in creating linkages
#  ## between fragments 
#
#  set topo([list core-repeat bapad]) {
#   {C2 C1 N1 C1 C2 C5} 
#   {C1 C2 N2 C1 C2 C5}
#  }
#
#  set topo([list repeat-repeat bapad]) {
#   {C2 C3 N1 C1 C2 C5}
#   {C2 C4 N2 C1 C2 C5}
#  }
#
#  set topo([list repeat-term
#
#}


proc ::DBTools::topo_bapad {} {

  variable topo

  ## Whitespace is ignored

  lappend topo bapad {

        name {
          core   EDH
          repeat AMP
          term   AHH 
        }

        fname { 
          core   EDH.mol2
          repeat AMP.mol2
          term   AHH.mol2
        }

        link {
          
          core-repeat {
            {C2 C1 N1 C1 C2 C5} 
            {C1 C2 N2 C1 C2 C5}
          }
          
          repeat-repeat {
            {C2 C3 N1 C1 C2 C5}
            {C2 C4 N2 C1 C2 C5}
          }  
        
          repeat-term {}
      
        }

        geometry {
          core-repeat {
            {175.0 1.35 125.7 -170 115.3 -76.0}
            {175.0 1.35 125.7 -170 115.3 -76.0}
          }

          repeat-repeat {
            {175.0 1.35 125.7 -170 115.3 -76.0}
            {175.0 1.35 125.7 -170 115.3 -76.0}
          }

          repeat-term {

          }
        }

        bonds {
          core-repeat {
            {N1 C1}
            {N2 C1}
          }

          repeat-repeat {
            {N1 C1}
            {N1 C1}
          }

          repeat-term {}

        }
  }

  return 

}

::DBTools::topo_bapad
