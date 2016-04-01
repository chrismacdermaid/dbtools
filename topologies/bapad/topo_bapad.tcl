if {![namespace exists ::DBTools]} {
  namespace eval ::DBTools:: {}
  set testflag 1
}

## {95.0 1.35 125.7 -170 115.3 -130.0}
## {95.5 1.35 125.7 -170 115.3 -130.0}

proc ::DBTools::topo_bapad {} {

  variable topo

  set resn bapad

  catch {
    if {[dict exists $::DBTools::topo $resn]} {
      dict unset $::DBTools::topo $resn
    }
  } 

  ## Whitespace is ignored

  lappend topo $resn {

        name {
          core   EDH
          repeat AMP
          term   AHH 
        }

        fname { 
          core   edh.mol2
          repeat repeat.pdb 
          term   ahh.mol2
        }

        link {
          
          core-repeat {
            {C2 C1 N1 C1 C2 O1} 
            {C1 C2 N2 C1 C2 O1}
          }
          
          repeat-repeat {
            {C2 C3 N1 C1 O1 C2}
            {C2 C4 N2 C1 O1 C2}
          }  
        
          repeat-term {
            {C2 C3 N1 C1 O1 C2}
            {C2 C4 N2 C1 O1 C2}
          }

          core-term {
            {C2 C1 N1 C1 O1 C2} 
            {C1 C2 N2 C1 O1 C2}
          }
        }

        geometry {
          core-repeat {
            {175.0 1.35 125.7 -170 115.3 178}
            {175.0 1.35 125.7 -170 115.3 178}
          }

          repeat-repeat {
            {175.0 1.35 125.7 0.0 120.0 180.0}
            {175.0 1.35 125.7 0.0 120.0 180.0}
          }

          repeat-term {
            {175.0 1.35 125.7 0.0 120.0 178.0}
            {175.0 1.35 125.7 0.0 120.0 178.0}
          }

          core-term {
            {175.0 1.35 125.7 0.0 120.0 178.0}
            {175.0 1.35 125.7 0.0 120.0 178.0}
          }
        }

        bonds {
          core-repeat {
            {N1 C1}
            {N2 C1}
          }

          repeat-repeat {
            {N1 C1}
            {N2 C1}
          }

          repeat-term {
            {N1 C1}
            {N2 C1}
          }

          core-term {
            {N1 C1}
            {N2 C1}
          }

        }
  }

  return 

}

::DBTools::topo_bapad
