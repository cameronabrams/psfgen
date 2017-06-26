# loopmc.tcl
# This set of TcL procedures are intended for use by VMD in a psfgen session 
# in which residues are modelled in using "residue" commands in the 
# segment stanzas.  The main procedure is "do_loop_mc" which
# uses a monte-carlo calculation to close loops.
#
# cameron f abrams, drexel u., 2017
# cfa22@drexel.edu
#
#
# generates a sequence of integers; use as [.. $i $j]
proc .. {from to} {
    set o {}
    if {$from <= $to} {
        for {set i $from} {$i <= $to} {incr i}    {lappend o $i}
    } else {
        for {set i $from} {$i >= $to} {incr i -1} {lappend o $i}
    }
    return $o
}

# a VERY rough energy calculation -- just counts contacts 
proc roughenergy {sel1 sel2 cut}  {
  set E 0.0
  set efac 1.0
  if { [$sel1 num] > 0 && [$sel2 num] > 0 } {
   set dat [measure contacts $cut $sel1 $sel2]
   set nc [llength [lindex $dat 0]]
   set E [expr $nc * $efac]
  }
  return $E
}

proc checknum { num msg } {
  if { $num == 0 } { 
    puts "$msg"
    #exit
  }
}

# rotates all atoms in chain c-terminal to residue r up to and 
# including residue rend in chain c around residue r's phi angle 
# by deg degrees in molecule with id molid.  Recall phi:  (-C)-N--CA-C
# Note:  the key in this procedure is the definition of the atomselection
# containing the atoms that will rotate.  We do NOT include
# in this set atoms N or HN in resid "r" nor atoms "C" and "O" in 
# resid "rend".  The former is because they are N-terminal to the
# phi bond, and the latter is because psfgen seems to build 
# the C and O atom positions based on the positions of N, CA, and CB
# in the NEXT residue, rather than on the internal coordinates
# of the built residue.  That means, when we use psfgen to build
# the raw loop from residue i to j inclusive, it is the CA-C bond
# residue j that is non-natural in length.
proc Crot_phi { r rend c molid deg } {
   set rot [atomselect $molid "((residue $r and not name N and not name HN) or (residue > $r and residue < $rend) or (residue $rend and not name C and not name O)) and chain $c"]; checknum [$rot num] "no rot in Crot_phi";
   set n [atomselect $molid "residue $r and name N"] ; checknum [$n num] "No N in Crot_phi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_phi";
   set p1 [lindex [$n get {x y z}] 0]
   set p2 [lindex [$ca get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $n delete
}

proc Crot_psi { r rend c molid deg } {
   set rot [atomselect $molid "((residue $r and name O) or (residue > $r and residue < $rend) or (residue $rend and not name C and not name O)) and chain $c"]; checknum [$rot num] "no rot in Crot_psi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_psi";
   set co [atomselect $molid "residue $r and name C"] ; checknum [$co num] "No C in Crot_psi";
   set p1 [lindex [$ca get {x y z}] 0]
   set p2 [lindex [$co get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $co delete
}

proc Crot_psi_special { r rend c molid deg } {
   set rot [atomselect $molid "((residue > $r and residue < $rend) or (residue $rend and not name C and not name O)) and chain $c"]; checknum [$rot num] "no rot in Crot_psi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_psi";
   set co [atomselect $molid "residue $r and name C"] ; checknum [$co num] "No C in Crot_psi";
   set p1 [lindex [$ca get {x y z}] 0]
   set p2 [lindex [$co get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $co delete
}

# rotates all residues C-terminal to "r" up to and including "rend" 
# (EXCEPT NOT C or O in rend) around the peptide bond between "r" and
# "r"+1.
proc Crot_omega { r rend c molid deg } {
   set rot [atomselect $molid "((residue > $r and residue < $rend) or (residue $rend and not name C and not name O)) and chain $c"]; checknum [$rot num] "no rot in Crot_omega";
   set co [atomselect $molid "residue $r and name C"] ; checknum [$co num] "No C in Crot_omega";
   set n [atomselect $molid "residue [expr $r + 1] and name N"] ; checknum [$n num] "No N in Crot_omega";
   set p1 [lindex [$co get {x y z}] 0]
   set p2 [lindex [$n get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $n delete
   $co delete
}

# returns a random integer between $min and $max, inclusive
proc irand_dom { min max } { 
   return [expr int(rand() * ($max - $min + 1)) + $min]
}

# uses a Metropolis MC algorithm w/ "temperature" to shrink the 
# CA-C bond in the last residue indicatedin the "residueList"
# by changing phi,psi angles in residues in residueList 
# on chain "c" in molecule molid.
# "k" is a spring constant for a fictitious spring 
# between the CA and C, and "r0" is the equilibrium length
# we'd like to establish.  "env" is the atomselection
# of the environment atoms which the loop atoms 
# are not allowed to overlap, where the overlap distance
# is "rcut" (A).
proc do_loop_mc { residueList c molid k r0 env rcut maxcycles temperature iseed } {

  set msel [atomselect $molid "chain $c and residue $residueList"]
  set mselnoh [atomselect $molid "chain $c and residue $residueList and noh"]
  set rend [lindex $residueList end]
  set nres [llength $residueList]

  set co [atomselect $molid "chain $c and residue $rend and name C"]
  set ca [atomselect $molid "chain $c and residue $rend and name CA"]
  set idx [list [$ca get index] [$co get index]]
  puts "$rend $idx"
  puts "CFALOOPMC) idx $idx distance C-CA($rend) [measure bond $idx]"

  expr srand($iseed)

  set nacc 0

  set SE [expr 0.5*$k*pow([measure bond $idx]-$r0,2)]
  set EE [roughenergy $mselnoh $env 3.0]
  set E [expr $SE + $EE]
  set E0 $E

  for {set cyc 0} { $cyc < $maxcycles } { incr cyc } {
    # save coordinates
    set SAVEPOS [$msel get {x y z}]
    for {set r 0} {$r < $nres} {incr r} {
      set i [irand_dom 0 [expr $nres-2]]
      set at [irand_dom 0 1]
      set av [expr 60 * [irand_dom 1 5]]
#      puts "CFALOOPMC) res $i angle-type $at angle-val $av"
      if { $at == 0 } { ; #phi
#        puts "CFALOOPMC) issuing Crot_phi $i $rend $c $molid $av"
        Crot_phi [lindex $residueList $i] $rend $c $molid $av
      } else {
#        puts "CFALOOPMC) issuing Crot_psi $i $rend $c $molid $av"
        Crot_psi [lindex $residueList $i] $rend $c $molid $av
      }
    }
    set SE [expr 0.5*$k*pow([measure bond $idx]-$r0,2)]
    set EE [roughenergy $mselnoh $env $rcut]
    set E [expr $SE + $EE]
    # decide to accept or reject this new conformation using a 
    # metropolis criterion
    set X [expr rand()]
    set arg [expr {($E0-$E)/$temperature}]
    if {[expr $arg < -20]} {
      set B 0.0
    } elseif {[expr $arg > 2.8]} {
      set B 1.1
    } else {
      # compute a Boltzmann factor
      set B [expr {exp($arg)}]
    }
    if {[expr {$X > $B}]} {
      # reject the move
      $msel set {x y z} $SAVEPOS
    } else {
      # accept the move
      set E0 $E
      incr nacc
      puts "CFALOOPMC) cycle $cyc nacc $nacc (ratio [format "%.5f" [expr (1.0*$nacc)/($cyc+1)]]) CA-C($rend): [format "%.2f" [measure bond $idx]] [format "link-penalty %.2f steric-penalty %.2f" $SE $EE]"
    }
  }
  puts "CFALOOPMC) cycle $cyc nacc $nacc (ratio [format "%.5f" [expr (1.0*$nacc)/($cyc+1)]]) CA-C($rend): [format "%.2f" [measure bond $idx]]"
}

# matrix inversion routine from http://wiki.tcl.tk/14921 (Keith Vetter)
proc Inverse3 {matrix} {
  if {[llength $matrix] != 3 ||
    [llength [lindex $matrix 0]] != 3 || 
    [llength [lindex $matrix 1]] != 3 || 
    [llength [lindex $matrix 2]] != 3} {
     error "wrong sized matrix"
  }
  set inv {{? ? ?} {? ? ?} {? ? ?}}
 
  # Get adjoint matrix : transpose of cofactor matrix
  for {set i 0} {$i < 3} {incr i} {
    for {set j 0} {$j < 3} {incr j} {
      lset inv $i $j [_Cofactor3 $matrix $i $j]
    }
  }
  # Now divide by the determinant
  set det [expr {double([lindex $matrix 0 0]   * [lindex $inv 0 0]
                 + [lindex $matrix 0 1] * [lindex $inv 1 0]
                 + [lindex $matrix 0 2] * [lindex $inv 2 0])}]
  if {$det == 0} {
    error "non-invertable matrix"
  }
    
  for {set i 0} {$i < 3} {incr i} {
    for {set j 0} {$j < 3} {incr j} {
      lset inv $i $j [expr {[lindex $inv $i $j] / $det}]
    }
  }
  return $inv
}

proc _Cofactor3 {matrix i j} {
  array set COLS {0 {1 2} 1 {0 2} 2 {0 1}}
  foreach {row1 row2} $COLS($j) break
  foreach {col1 col2} $COLS($i) break
    
  set a [lindex $matrix $row1 $col1]
  set b [lindex $matrix $row1 $col2]
  set c [lindex $matrix $row2 $col1]
  set d [lindex $matrix $row2 $col2]
 
  set det [expr {$a*$d - $b*$c}]
  if {($i+$j) & 1} { set det [expr {-$det}]}
  return $det
}

# returns the position of an amide N given the upstream position
# of a CA, C, and O.  This is necessary because psfgen does not know
# how to use IC's to generate de novo the position of the N of
# a modeled in residue off of a residue with an arbitrary position.
proc cacoIn_nOut { resid chain molid } {
  set rn {? ? ?}

  set ca [atomselect $molid "chain $chain and resid $resid and name CA"]
  set c  [atomselect $molid "chain $chain and resid $resid and name C"]
  set o  [atomselect $molid "chain $chain and resid $resid and name O"]

  set r1 [lindex [$ca get {x y z}] 0]
  set r2 [lindex [$c  get {x y z}] 0]
  set r3 [lindex [$o  get {x y z}] 0]

  set R21 [vecsub $r2 $r1]
  set r21 [vecscale $R21 [expr 1.0/[veclength $R21]]]
  set R32 [vecsub $r3 $r2]
  set r32 [vecscale $R32 [expr 1.0/[veclength $R32]]]
  set c [veccross $r21 $r32]

  set MAT { {? ? ?} {? ? ?} {? ? ?} }
  lset MAT 0 0 [lindex $c 0]
  lset MAT 0 1 [lindex $c 1]
  lset MAT 0 2 [lindex $c 2]
  lset MAT 1 0 [lindex $r21 0]
  lset MAT 1 1 [lindex $r21 1]
  lset MAT 1 2 [lindex $r21 2]
  lset MAT 2 0 [lindex $r32 0]
  lset MAT 2 1 [lindex $r32 1]
  lset MAT 2 2 [lindex $r32 2]

#  puts "MAT $MAT"
#  exit
  set b { ? ? ? }
  lset b 0 0
  lset b 1 [expr -cos(3.14159*114.44/180.0)] 
  lset b 2 [expr cos(3.14159*123.04/180.0)]
  
  set AMAT [Inverse3 $MAT]

  set rnd { ? ? ? }
  lset rnd 0 [expr [lindex $AMAT 0 0] * [lindex $b 0] + [lindex $AMAT 0 1] * [lindex $b 1] + [lindex $AMAT 0 2] * [lindex $b 2]]
  lset rnd 1 [expr [lindex $AMAT 1 0] * [lindex $b 0] + [lindex $AMAT 1 1] * [lindex $b 1] + [lindex $AMAT 1 2] * [lindex $b 2]]
  lset rnd 2 [expr [lindex $AMAT 2 0] * [lindex $b 0] + [lindex $AMAT 2 1] * [lindex $b 1] + [lindex $AMAT 2 2] * [lindex $b 2]]

  set rn [vecadd $r2 [vecscale $rnd 1.355]]

  return $rn
}
