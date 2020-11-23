# The purpose of this script is to define the 'axeq' procedure.  The axeq procedure
# assists in determining whether the bond to a given glycan's C1 atom is 'equatorial'
# or 'axial'.  The operations to determine this are based solely on the geometry of 
# the glycan ring containing the C1 and the position of the atom from the 'incoming'
# ring that connects to it.
#
# Cameron Abrams cfa22@drexel.edu

proc axeq { ose_resid molid chain in_name c1_resid } {
   set resname [[atomselect $molid "resid $ose_resid and chain $chain and name C1"] get resname]
   if { $c1_resid != -1 } {
      set nra_sel [atomselect $molid "(resid $ose_resid and chain $chain and not name C1 C2 C3 C4 C5 O5) or (resid $c1_resid and chain $chain and not name C1 C2 C3 C4 C5 O5)"]
   } else {
      set nra_sel [atomselect $molid "resid $ose_resid and chain $chain and not name C1 C2 C3 C4 C5 O5"]
   }
   set nra_i [$nra_sel get index]
   set nra_n [$nra_sel get name]
   set nra_x [$nra_sel get x]
   set nra_y [$nra_sel get y]
   set nra_z [$nra_sel get z]
   set ra_sel [atomselect $molid "resid $ose_resid and chain $chain and name C1 C2 C3 C4 C5 O5"]
   set ra_i [$ra_sel get index]
   set ra_n [$ra_sel get name]
   set ra_x [$ra_sel get x]
   set ra_y [$ra_sel get y]
   set ra_z [$ra_sel get z]

   foreach n $nra_i nn $nra_n nx $nra_x ny $nra_y nz $nra_z {
      set pos($nn) [list $nx $ny $nz]
   }
   set i 0
   foreach r $ra_i rn $ra_n rx $ra_x ry $ra_y rz $ra_z {
      set pos($rn) [list $rx $ry $rz]
      set forp($rn) [lindex $ra_n [expr ($i+2)%6]]
      set bakp($rn) [lindex $ra_n [expr ($i-2)%6]]
      foreach n $nra_i nn $nra_n nx $nra_x ny $nra_y nz $nra_z {
	      set tb [measure bond [list $r $n]]
	      if {$tb < 2.0} {
		      set ligand_by_name($rn) $nn
		      set ligand_by_index($r) $n
         }
      }
      incr i
   }
   foreach {rn ln} [array get ligand_by_name] {
	   set forvec [vecsub $pos($rn) $pos($forp($rn))]
	   set bakvec [vecsub $pos($rn) $pos($bakp($rn))]
	   set pcross [veccross $forvec $bakvec]
	   set ligvec [vecsub $pos($rn) $pos($ln)]
	   set ligpdot [expr abs([vecdot $ligvec $pcross])]
	   if { $ligpdot > 5.5 } {
	      set ligand_axeq($rn) "a"
	      set ligand_axeq($ln) "a"
	   } else {
         set ligand_axeq($rn) "b"
	      set liband_axeq($ln) "b"
	   }
      puts "ring atom $rn has forp $forp($rn) and bakp $bakp($rn) and ligand $ln bondlength [veclength $ligvec] ligpdot $ligpdot axeq $ligand_axeq($rn)"
   }
   return $ligand_axeq($in_name)
}
