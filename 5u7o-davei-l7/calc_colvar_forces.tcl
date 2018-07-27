proc my_avoid_force { d } {
   if { [expr $d < 20.0 ] } {
      return [expr -100.0 * ( $d - 20 ) ]
   } else {
      return 0.0
   }
}
proc calc_colvar_forces { ts } {
   set d1 [cv colvar trp3_1av value]
   set d2 [cv colvar trp3_2av value]
   set d3 [cv colvar trp3_3av value]
   set f1 [my_avoid_force $d1]
   set f2 [my_avoid_force $d2]
   set f3 [my_avoid_force $d3]
   if { [expr int(fmod($ts,100))] == "0" } { 
     puts "$ts [format "%.5f" $d1] [format "%.5f" $d2] [format "%.5f" $d3] [format "%.5f" $f1] [format "%.5f" $f2] [format "%.5f" $f3]"
   }
   cv colvar trp3_1av addforce $f1
   cv colvar trp3_2av addforce $f2
   cv colvar trp3_3av addforce $f3]
}

