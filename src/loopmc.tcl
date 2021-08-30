# loopmc.tcl
# This set of TcL procedures are intended for use by VMD in a psfgen session 
# in which residues are modelled in using "residue" commands in the 
# segment stanzas.  The main procedure is "do_loop_mc" which
# uses a monte-carlo calculation to close loops.
#
# cameron f abrams, drexel u., 2017-2020
# cfa22@drexel.edu
#
# some C functions that help with quick bond rotations
if {![info exists PSFGEN_BASEDIR]} {
    if {[info exists env(PSFGEN_BASEDIR)]} {
        set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
    } else {
        set PSFGEN_BASEDIR $env(HOME)/research/psfgen
    }
}
source ${PSFGEN_BASEDIR}/src/cfarot.tcl

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

proc pc3 { pt } {
  return "( [lindex $pt 0] , [lindex $pt 1] , [lindex $pt 2] )"
}

proc my_increment { numlet } {
   if { [string is integer $numlet] == "1" } {
     return [expr $numlet + 1]
   }
   set c [format %c [expr [scan [string index $numlet end] %c] + 1]]
   set nn [string replace $numlet end end $c]
   return $nn
}

# computes overlap energy between atoms in sel1 and sel2.  "sel2" is assumed to be a static background that
# does not include the atoms in sel1.  Atom indices in sel1 are assumed not to change.
proc roughenergy_setup { sel1 sel2 cellsize } {
  global _r1
  global _r2
  global _x1
  global _y1
  global _z1
  global _x2
  global _y2
  global _z2
  global _n1
  global _n2
  set _r1 [intListToArray [$sel1 get index]]
  set _r2 [intListToArray [$sel2 get index]]
  set _x1 [ListToArray [$sel1 get x]]
  set _y1 [ListToArray [$sel1 get y]]
  set _z1 [ListToArray [$sel1 get z]]
  set _x2 [ListToArray [$sel2 get x]]
  set _y2 [ListToArray [$sel2 get y]]
  set _z2 [ListToArray [$sel2 get z]]
  set _n1 [$sel1 num]
  set _n2 [$sel2 num]
  #puts "Calling my_roughenergy_setup"
  set ls [my_roughenergy_setup $_x2 $_y2 $_z2 $_n2 $cellsize]
  return $ls
}

proc roughenergy_cleanup { ls } {
  my_roughenergy_cleanup $ls
}

proc roughenergy { sel1 cut sigma epsilon shift bs ls }  {
  global _r1
  global _r2
  global _x1
  global _y1
  global _z1
  global _x2
  global _y2
  global _z2
  global _n1
  global _n2

  set E 0.0

  if { $_n1 == 0 || $_n2 == 0 } { return $E }
  # update positions in 1
  ListToArray_Data $_x1 [$sel1 get x]
  ListToArray_Data $_y1 [$sel1 get y]
  ListToArray_Data $_z1 [$sel1 get z]
  set E [my_roughenergy $_r1 $_x1 $_y1 $_z1 $_n1 $_r2 $_n2 $cut $sigma $epsilon $shift $bs $ls]
  
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

# same as above, but treats last residue as C-terminus
proc Crot_phi_toCterm { r rend c molid deg } {
   set rot [atomselect $molid "((residue $r and not name N and not name HN) or (residue > $r and residue <= $rend)) and chain $c"]; checknum [$rot num] "no rot in Crot_phi";
   set n [atomselect $molid "residue $r and name N"] ; checknum [$n num] "No N in Crot_phi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_phi";
   set p1 [lindex [$n get {x y z}] 0]
   set p2 [lindex [$ca get {x y z}] 0]
   set ax [vecsub $p2 $p1]
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

proc Crot_psi_toCterm { r rend c molid deg } {
   set rot [atomselect $molid "((residue $r and name O) or (residue > $r and residue <= $rend)) and chain $c"]; checknum [$rot num] "no rot in Crot_psi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_psi";
   set co [atomselect $molid "residue $r and name C"] ; checknum [$co num] "No C in Crot_psi";
   set p1 [lindex [$ca get {x y z}] 0]
   set p2 [lindex [$co get {x y z}] 0]
   set ax [vecsub $p2 $p1]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $co delete
}

proc range { from to step } {
  set result {}
  for {set i $from} {$i<=$to} {incr i $step} {
    lappend result $i
  }
  return $result
}

proc lay_loop { molid c loop maxcycles } {
  set nr [llength $loop]
  set loopsel [atomselect $molid "chain $c and resid [join $loop]"]
  set residue_numbers [[atomselect $molid "[$loopsel text] and name CA"] get residue]
  set env [atomselect $molid "same residue as exwithin 4.0 of (chain $c and resid [join $loop])"]
  set residuenum_end [lindex $residue_numbers end]
  for { set i 0 } { $i < $nr } { incr i } {
    # rotate phi angle and psi angle to minimize number of contacts between residue and 
    # its environment
    set rsel [atomselect $molid "chain $c and resid [lindex $loop $i] to [lindex $loop end]"]
    set residuenum1 [lindex $residue_numbers $i]
    set CON_STRUCT [measure contacts 2.0 $rsel $env]
    set CON [llength [lindex $CON_STRUCT 0]]
#    puts "LAYLOOP) ${c}[lindex $loop $i] INIT $CON"
    for { set t 0 } { $t < $maxcycles && $CON > 0 } { incr t } {
      set SAVEPOS [$loopsel get {x y z}]
      set rphi [expr (1-2*rand())*120.0]
      #set rpsi [expr (1-2*rand())*60.0]
      Crot_phi_toCterm $residuenum1 $residuenum_end $c $molid $rphi
      if { $i > [expr $nr - 1] } {
        Crot_psi_toCterm $residuenum1 $residuenum_end $c $molid $rphi
      }
      $env update
      set TRICON_STRUCT [measure contacts 2.0 $rsel $env]
      set TRICON  [llength [lindex $TRICON_STRUCT 0]]
      if { [expr $TRICON < $CON] } {
        # accept this move
        set CON $TRICON
        puts "LAYLOOP) ${c}:[lindex $loop 0]-[lindex $loop $i] $t $CON"
      } else {
        # reject this move
        $loopsel set {x y z} $SAVEPOS
      }
    }
    $rsel delete
  }
  $env delete
  $loopsel delete
}

# rotate the side chain of residue r of chain c in mol molid around
# chi1 by deg degrees
proc SCrot_chi1 { r c molid deg } {
   set rot [atomselect $molid "residue $r and not name N HN CA CB HA C O"]; checknum [$rot num] "no rot in SCrot_chi1";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "no CA in SCrot_chi1";
   set cb [atomselect $molid "residue $r and name CB"]; checknum [$cb num] "no CB in SCrot_chi1";
   set p1 [lindex [$ca get {x y z}] 0]
   set p2 [lindex [$cb get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $cb delete
}

# rotate the side chain of residue r of chain c in mol molid around
# chi2 by deg degrees
proc SCrot_chi2 { r c molid deg } {
   set rot [atomselect $molid "residue $r and not name N HN CA CB HA C O HB1 HB2"]; checknum [$rot num] "no rot in SCrot_chi1";
   set cb [atomselect $molid "residue $r and name CB"]; checknum [$cb num] "no CB in SCrot_chi1";
   set cg [atomselect $molid "residue $r and name CG"]; checknum [$cg num] "no CG in SCrot_chi1";
   set p1 [lindex [$cb get {x y z}] 0]
   set p2 [lindex [$cg get] {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $cb delete
   $cg delete
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

# logs frames to a logging molecule
proc log_addframe { molid logid } {
   if { $logid != "-1" } {
     animate dup $logid
     [atomselect $logid all] set x [[atomselect $molid all] get x]
     [atomselect $logid all] set y [[atomselect $molid all] get y]
     [atomselect $logid all] set z [[atomselect $molid all] get z]
#     [atomselect $molid all] writepdb "tmp.pdb"
#     animate read pdb tmp.pdb $logid
     puts "Molid $molid - logging molecule $logid has [molinfo $logid get numframes] frames."
#     exec rm -f tmp.pdb
   }
}

# uses a Metropolis MC algorithm w/ "temperature" to shrink the 
# CA-C bond in the last residue indicated in the "residueList"
# by changing phi,psi angles in residues in residueList 
# on chain "c" in molecule molid.
# "k" is a spring constant for a fictitious spring 
# between the CA and C, and "r0" is the equilibrium length
# we'd like to establish.  "env" is the atomselection
# of the environment atoms which the loop atoms 
# are not allowed to overlap, where the overlap distance
# is "rcut" (A).
# "logid" is the optional molecule id of a logging molecule (-1 means do nothing)
proc do_loop_mc { residueList c molid k r0 env sigma epsilon rcut maxcycles temperature iseed logid logevery } {

  set msel [atomselect $molid "chain $c and residue $residueList"]
  set mselnoh [atomselect $molid "chain $c and residue $residueList and noh"]
  set exind [$mselnoh get index]
  set envex [atomselect $molid "[$env text] and not index $exind"]

  if { [$envex num] == 0 } {
    puts "Error: environment minus loop has zero atoms."
    return
  }

  set bs [make_bondstruct $molid $mselnoh]
  #bondstruct_print $bs
#  set mselca [atomselect $molid "chain $c and residue $residueList and name CA"]
#  set envca [atomselect $molid "([$env text]) and name CA"]

  set rend [lindex $residueList end]
  set nres [llength $residueList]

  set co [atomselect $molid "chain $c and residue $rend and name C"]
  set ca [atomselect $molid "chain $c and residue $rend and name CA"]
  set idx [list [$ca get index] [$co get index]]
  set ibl [format "%.2f" [measure bond $idx]]
  puts "CFALOOPMC) residue $rend CA_[$ca get index]-C_[$co get index]: $ibl A"

  expr srand($iseed)

  set nacc 0
  
  #puts "roughenergy setup..."
  set ls [roughenergy_setup $mselnoh $envex $rcut]
  #puts "calc ($mselnoh) ($rcut) ($sigma) ($epsilon) ($bs) ($ls)..."
  set SE [expr 0.5*$k*pow([measure bond $idx]-$r0,2)]
  
  set EE [roughenergy $mselnoh $rcut $sigma $epsilon $bs $ls]
  set E [expr $SE + $EE]
  set E0 $E
  #puts "EE $EE"
  flush stdout

  for {set cyc 0} { $cyc < $maxcycles } { incr cyc } {
    # save coordinates
    set SAVEPOS [$msel get {x y z}]
    # do a cycle of $nres rotation attempts
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
    #set EE [roughenergy $mselnoh $env $rcut]
    set EE [roughenergy $mselnoh $rcut $sigma $epsilon $bs $ls]
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
      puts "CFALOOPMC) ($rend) cyc $cyc na $nacc ([format "%.5f" [expr (1.0*$nacc)/($cyc+1)]]) CA-C: [format "%.2f" [measure bond $idx]] [format "lnk-pnlty %.2f strc-pnlty %.2f" $SE $EE]"
      if { [expr $nacc % $logevery == 0 ] } {
        log_addframe $molid $logid
      }
    }
  }
  puts "CFALOOPMC) ($rend) cyc $cyc na $nacc ([format "%.5f" [expr (1.0*$nacc)/($cyc+1)]]) CA-C: [format "%.2f" [measure bond $idx]]"
  free_bondstruct $bs
  roughenergy_cleanup $ls
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

#  puts "$molid $chain $resid [$ca get name] [$c get name] [$o get name]"

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

# returns the 3-vector cartesian coordinates of the center of the 
# circle that passes through these three points.
proc center_of { pt1 pt2 pt3 } {
  set ax [lindex $pt1 0]
  set ay [lindex $pt1 1]
  set az [lindex $pt1 2]
  set bx [lindex $pt2 0]
  set by [lindex $pt2 1]
  set bz [lindex $pt2 2]
  set cx [lindex $pt3 0]
  set cy [lindex $pt3 1]
  set cz [lindex $pt3 2]

  set yzzy [expr ($ay-$by)*($bz-$cz) - ($az-$bz)*($by-$cy)]
  set zxxz [expr ($az-$bz)*($bx-$cx) - ($ax-$bx)*($bz-$cz)]
  set xyyx [expr ($ax-$bx)*($by-$cy) - ($ay-$by)*($bx-$cx)]

  set ax2 [expr $ax * $ax]
  set ay2 [expr $ay * $ay]
  set az2 [expr $az * $az]
  set bx2 [expr $bx * $bx]
  set by2 [expr $by * $by]
  set bz2 [expr $bz * $bz]
  set cx2 [expr $cx * $cx]
  set cy2 [expr $cy * $cy]
  set cz2 [expr $cz * $cz]

  set MAT { {? ? ?} {? ? ?} {? ? ?} }
  lset MAT 0 0 [expr 2*($bx-$ax)]
  lset MAT 0 1 [expr 2*($by-$ay)]
  lset MAT 0 2 [expr 2*($bz-$az)]
  lset MAT 1 0 [expr 2*($cx-$bx)]
  lset MAT 1 1 [expr 2*($cy-$by)]
  lset MAT 1 2 [expr 2*($cz-$bz)]
  lset MAT 2 0 $yzzy
  lset MAT 2 1 $zxxz
  lset MAT 2 2 $xyyx
 
  set BB { ? ? ? }
  lset BB 0 [expr $bx2-$ax2+$by2-$ay2+$bz2-$az2]
  lset BB 1 [expr $cx2-$bx2+$cy2-$by2+$cz2-$bz2]
  lset BB 2 [expr $ax*$yzzy+$ay*$zxxz+$az*$xyyx]

  set AMAT [Inverse3 $MAT]

  set o { ? ? ? }
  lset o 0 [expr [lindex $AMAT 0 0] * [lindex $BB 0] + [lindex $AMAT 0 1] * [lindex $BB 1] + [lindex $AMAT 0 2] * [lindex $BB 2]]
  lset o 1 [expr [lindex $AMAT 1 0] * [lindex $BB 0] + [lindex $AMAT 1 1] * [lindex $BB 1] + [lindex $AMAT 1 2] * [lindex $BB 2]]
  lset o 2 [expr [lindex $AMAT 2 0] * [lindex $BB 0] + [lindex $AMAT 2 1] * [lindex $BB 1] + [lindex $AMAT 2 2] * [lindex $BB 2]]

  #puts "check: [veclength [vecsub $pt1 $o]] [veclength [vecsub $pt2 $o]] [veclength [vecsub $pt3 $o]]"

  return $o

}

# generic molecular fragment rotation
# "sel" is an existing atomselection for the molecule
# i and j are the indices of the two bonded atoms
# that define the rotation.  It is assumed that
# i and j share a rotatable bond; that is,
# if the i-j bond were broken, the molecule would be
# broken into two noncontiguous fragments.  The
# rotation is about the i-j bond and includes all
# atoms on the "j"-side of the i-j bond.
# cfa 2018
proc genbondrot { molid sel i j deg } {
   set ilist [$sel get index]
   set blist [$sel getbonds]
   set ii [lsearch $ilist $i]
   set jj [lsearch $ilist $j]
   set xlist [$sel get x]
   set ylist [$sel get y]
   set zlist [$sel get z]
   set ix [lindex $xlist $ii]
   set iy [lindex $ylist $ii]
   set iz [lindex $zlist $ii]
   set jx [lindex $xlist $jj]
   set jy [lindex $ylist $jj]
   set jz [lindex $zlist $jj]
   set bi [lindex $blist $ii]
   set isj [lsearch $bi $j]
   if { $isj == -1 } { 
       puts "ERROR: $j not found in $i's bondlist $bi"
   } else {
     set rlist [list $j]
     set bj [lindex $blist $jj]
     foreach n $bj {
        if { $n != $i } {
           lappend rlist $n
        }
     }
     set grow 1
     while { $grow } {
       set grow 0
       set rfrag {}
       foreach k $rlist {
          set kk [lsearch $ilist $k]
          set bk [lindex $blist $kk]
          foreach kn $bk {
             if { $kn != $i } {
               set pr [expr [lsearch $rlist $kn] + [lsearch $rfrag $kn]]
               if { $pr == -2 } {
                 lappend rfrag $kn
               }
             } 
          }
       }
       if { [llength $rfrag] > 0 } {
         set rlist [concat $rlist $rfrag]
         set grow 1
       }
     }
     set rsel [atomselect $molid "index $rlist"]
     $rsel move [trans bond [list $ix $iy $iz] [list $jx $jy $jz] $deg degrees]
     $rsel delete
   }
}

proc fold_alpha_helix { molid sel ebfn } {
   set ni {}
   set ci {}
   set cai {}
#   puts "molid $molid"
   foreach an [$sel get name] i [$sel get index] {
      if { $an == "C" } {
         lappend ci $i
      } elseif { $an == "N" } {
         lappend ni $i
      } elseif { $an == "CA" } {
         lappend cai $i
      }
   }
   set nn [llength $ni]
   set nc [llength $ci]
   set nca [llength $cai]
#   puts "N $nn $ni"
#   puts "CA $nca $cai"
#   puts "C $nc $ci"
   set in 0
   set ica 0
   set ic 0
   if { $ebfn != "0" } {
     set fp [open $ebfn "w"]
   }
   while { $in < [expr $nn - 1] } {
     set psi {}
     set phi {}
     lappend psi [lindex $ni $in]
     lappend psi [lindex $cai $ica]
     lappend psi [lindex $ci $ic]
     lappend phi [lindex $ci $ic]
     set in [expr $in + 1]
     lappend psi [lindex $ni $in]
     lappend phi [lindex $ni $in]
     set ica [expr $ica + 1]
     lappend phi [lindex $cai $ica]
     set ic [expr $ic + 1]
     lappend phi [lindex $ci $ic]
     set PSIM [measure dihed $psi]
     set PHIM [measure dihed $phi]
#     puts "psi $psi ($PSIM) phi $phi ($PHIM)"
     set SPHI [expr -60 - ($PHIM)]
     set SPSI [expr -45 - ($PSIM)]
     if { $ebfn != "0" } {
       puts $fp "dihedral $psi 10.0 -45.0"
       puts $fp "dihedral $phi 10.0 -60.0"
     }
     genbondrot $molid $sel [lindex $psi 1] [lindex $psi 2] $SPSI
     genbondrot $molid $sel [lindex $phi 1] [lindex $phi 2] $SPHI
#     puts "$SPHI $SPSI"
   }
   if { $ebfn != "0" } {
     close $fp
   }
}

proc random_loop { molid sel } {
   set ni {}
   set ci {}
   set cai {}
#   puts "molid $molid"
   foreach an [$sel get name] i [$sel get index] {
      if { $an == "C" } {
         lappend ci $i
      } elseif { $an == "N" } {
         lappend ni $i
      } elseif { $an == "CA" } {
         lappend cai $i
      }
   }
   set nn [llength $ni]
   set nc [llength $ci]
   set nca [llength $cai]
#   puts "N $nn $ni"
#   puts "CA $nca $cai"
#   puts "C $nc $ci"
   set in 0
   set ica 0
   set ic 0
   while { $in < [expr $nn - 1] } {
     set psi {}
     set phi {}
     lappend psi [lindex $ni $in]
     lappend psi [lindex $cai $ica]
     lappend psi [lindex $ci $ic]
     lappend phi [lindex $ci $ic]
     set in [expr $in + 1]
     lappend psi [lindex $ni $in]
     lappend phi [lindex $ni $in]
     set ica [expr $ica + 1]
     lappend phi [lindex $cai $ica]
     set ic [expr $ic + 1]
     lappend phi [lindex $ci $ic]
     set PSIM [measure dihed $psi]
     set PHIM [measure dihed $phi]
     puts "psi $psi ($PSIM) phi $phi ($PHIM)"
     # pick a random (phi,psi) from ramachandran plot
     set PHIR [irand_dom -340 -60]
     set PSIR [irand_dom -60 300]
     set SPHI [expr ($PHIR) - ($PHIM)]
     set SPSI [expr ($PSIR) - ($PSIM)]
     genbondrot $molid $sel [lindex $psi 1] [lindex $psi 2] $SPSI
     genbondrot $molid $sel [lindex $phi 1] [lindex $phi 2] $SPHI
#     puts "$SPHI $SPSI"
   }
}

# do_flex_mc: loop MC on atoms in msel.  
# molid is the molcule id.
# msel is the selection of atoms that contain all rotatable bonds.
# ri is a list of atom indices in msel that are rotatable bond "left" partners.
# rj is a list of atom indices in msel that are rotatable bond "right" partners.
# ri and rj must be set by the calling routine!
# fa is the index of an atom that defines an end of the molecule that is not
#   permitted to move due to any bond rotation (if -1, all bonds rotate)
# All bonds in which a left partner is on ri and a right partner is on rj will
#   be rotated, with atoms on the right of the bond moving and the atoms on the left
#   of the bond fixed (this is the meaning of "left" and "right")
# "i" and "j" are reference atom indices in msel that act as an attractor; it is assumed
#   that i is on the loop and j is not; k is the spring constant for a harmonic attractor.
# envsel is a second atomselection of atoms that form the "environment" for which
#   overlaps are not allowed (typically all) and rcut is the minimum distance between
#   any pair of atoms between the msel and envsel that would constitute a contact.
# maxcycles is the max number of mc cycles, where one cycle is a move of all rotatable
#   bonds by random amounts. 
# temperature is the Metropolis temperature.
# iseed is the rng seed.
proc do_flex_mc { molid msel envsel refatominddict paramsdict iseed logid logevery logsaveevery } {

   upvar 1 $refatominddict refatoms
   upvar 1 $paramsdict params

   set mselnoh [atomselect $molid "([$msel text]) and noh"]
   #set bl [$msel getbonds]
   #set il [$mselnoh get index]
   set fa [dict get $refatoms fa]
   set i [dict get $refatoms ca]
   set j [dict get $refatoms c]
   set bs [make_bondstruct $molid $msel]
   bondstruct_deactivate_by_fixed $bs $fa
  # bondstruct_print $bs
   set exind [$msel get index]
   set envex [atomselect $molid "[$envsel text] and not index $exind"]
   puts "CFAFLEXMC) msel [$msel num] envex [$envex num] fa $fa"
   flush stdout
   if { $i != $j } { 
     puts "CFAFLEXMC) Initial attractor distance [format "%.2f" [measure bond [list $i $j]]] A"
   }
   set maxcycles [dict get $params nc]
   set dstop  [dict get $params dstop]
   set sstop  [dict get $params sstop]
   set k  [dict get $params mck]
   set temperature  [dict get $params temperature]
   set sigma  [dict get $params sigma]
   set epsilon [dict get $params epsilon]
   set rcut [dict get $params rcut]
   set maxanglestep [dict get $params maxanglestep]

   puts "CFAFLEXMC) Max cycles $maxcycles dattr-thresh $dstop strc-thresh $sstop k $k"
   puts "CFAFLEXMC) [bondstruct_getnrb $bs] rotatable bonds"
   puts "CFAFLEXMC) MC-Temperature $temperature sigma $sigma epsilon $epsilon"
   puts "CFAFLEXMC) Maximum angle displacement: $maxanglestep degrees"
   flush stdout

   set maxanglestep [expr $maxanglestep / 10.0]

   expr srand($iseed)

   set nacc 0

   set SE 0.0
   set dattr 0.0
   if { $i != $j } {
     set dattr [measure bond [list $i $j]]
     set SE [expr 0.5*$k*pow($dattr,2)]
   }
   set ls [roughenergy_setup $mselnoh $envex $rcut]
  #puts "calc ($mselnoh) ($rcut) ($sigma) ($epsilon) ($bs) ($ls)..."
   set EE [roughenergy $mselnoh $rcut $sigma $epsilon $bs $ls]
   set E [expr $SE + $EE]
   set lastEE $EE
   set lastSE $SE
   set E0 $E
   #puts "CFAFLEXMC) E0 $E0"
   set keep_cycling 1
   if { $EE < $sstop && $dattr < $dstop } {
      set keep_cycling 0
   }
   for {set cyc 0} { $cyc < $maxcycles && $keep_cycling == 1 } { incr cyc } {
      # save coordinates
      set SAVEPOS [$msel get {x y z}]
      set nrot 0
      for {set r 0} {$r < [bondstruct_getnrb $bs] } {incr r} {
         #set av [expr 60 * [irand_dom 1 5]]
         set av [expr $maxanglestep * [irand_dom -10 10]]
        # puts "cyc $cyc bond $r deg $av"
        set rr [bondstruct_r2b $bs $r]
         if { [bondstruct_isactive $bs $rr] } {
           bondrot_by_index $bs $molid $rr $av
           set nrot [expr $nrot + 1]
         }
      }
    
      if { $nrot == 0 } {
         puts "ERROR: no rotations performed"
         exit
      }
      if { $i != $j } {
        set dattr [measure bond [list $i $j]]
        set SE [expr 0.5*$k*pow($dattr,2)]
      } else {
        set SE 0.0
      }
      set EE [roughenergy $mselnoh $rcut $sigma $epsilon $bs $ls]
      set E [expr $SE + $EE]
     # puts " ... E $E"
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
        puts -nonewline "CFAFLEXMC) cyc $cyc na $nacc ([format "%.5f" [expr (1.0*$nacc)/($cyc+1)]]) "
        if { $i != $j } {
          puts -nonewline "attr dst: [format "%.2f" [measure bond [list $i $j]]] [format "attr-pnlty %.2f " $SE]"
        }
        puts "[format "strc-pnlty %.2f" $EE]"
        if { [expr $nacc % $logevery == 0 ] } {
          log_addframe $molid $logid
          if { [expr $nacc % $logsaveevery == 0] } {
             set loga [atomselect $logid all]
             animate write dcd "tmp.dcd" waitfor all sel $loga $logid
          }
        }
        set lastEE $EE
        set lastSE $SE
        if { $EE < $sstop && $dattr < $dstop } {
          set keep_cycling 0
        }
      }
   }
   puts -nonewline "CFAFLEXMC) Stop at cycle $cyc: "
   if { $i != $j } {
     puts -nonewline "attr dst: [format "%.2f" [measure bond [list $i $j]]] [format "attr-pnlty %.2f " $lastSE]"
   }
   puts "[format "strc-pnlty %.2f" $lastEE]"
   free_bondstruct $bs
   roughenergy_cleanup $ls
}

proc do_multiflex_mc { molid rotsel refatominddict paramsdict iseed logid logevery logsaveevery } {

   upvar 1 $refatominddict refatoms
   upvar 1 $paramsdict params

   set rotnoh [atomselect $molid "([$rotsel text]) and noh"]
   #set bl [$msel getbonds]
   #set il [$mselnoh get index]
   set falist [dict get $refatoms fa]
   set ilist [dict get $refatoms i]
   set jlist [dict get $refatoms j]
   set loc [lsearch $ilist -1]
   set trunc_ilist [lrange $ilist 0 [expr $loc-1]]
   set trunc_jlist [lrange $jlist 0 [expr $loc-1]]
 
   puts "CFAFLEXMC) Making bondstruct..."
   flush stdout 
   set bs [make_bondstruct $molid $rotsel $falist]

   # remove movable atoms from the background
   set exind [$rotsel get index]
   set envex [atomselect $molid "noh and not index $exind"]

   puts "CFAFLEXMC) rotsel [$rotsel num] envex [$envex num] falist $falist"
   set maxcycles [dict get $params nc]
   set dstop  [dict get $params dstop]
   set sstop  [dict get $params sstop]
   set k  [dict get $params mck]
   set temperature  [dict get $params temperature]
   set maxanglestep [dict get $params maxanglestep]

   set ljsigma  [dict get $params ljsigma]
   set ljepsilon [dict get $params ljepsilon]
   set ljcutoff  [dict get $params ljcutoff]
   set ljshift [dict get $params ljshift]
   set cellsize [dict get $params cellsize]

   puts "CFAFLEXMC) Max cycles $maxcycles dattr-thresh $dstop strc-thresh $sstop k $k"
   puts "CFAFLEXMC) Number of rotatable bonds [bondstruct_getnrb $bs]"
   puts "CFAFLEXMC) MC-Temperature $temperature"
   puts "CFAFLEXMC) sigma $ljsigma epsilon $ljepsilon cutoff $ljcutoff shift $ljshift"
   puts "CFAFLEXMC) Maximum angle displacement: $maxanglestep degrees"
   flush stdout

   puts "CFAFLEXMC) ilist $ilist"
   puts "CFAFLEXMC) jlist $jlist"
   flush stdout
   foreach i $ilist j $jlist {
      if { $i != $j } { 
        puts "CFAFLEXMC) Initial ($i)-($j) distance [format "%.2f" [measure bond [list $i $j]]] A"
        flush stdout
      }
   }

   set anglegradationfactor 10
   set maxanglestep [expr $maxanglestep / $anglegradationfactor]
   expr srand($iseed)
   set nacc 0
   set SE 0.0
   foreach i $ilist j $jlist {
      if { $i != $j } {
        set dattr [measure bond [list $i $j]]
        set SE [expr $SE+0.5*$k*pow($dattr,2)]
      }
   }
   set ls [roughenergy_setup $rotnoh $envex $cellsize]
  #puts "calc ($mselnoh) ($rcut) ($sigma) ($epsilon) ($bs) ($ls)..."
   puts -nonewline "roughenergy..."
   flush stdout
   set start [clock microseconds]
   set EE [roughenergy $rotnoh $ljcutoff $ljsigma $ljepsilon $ljshift $bs $ls]
   set stop [clock microseconds]
   puts "...done [expr ($stop-$start)/1.e6] s"
   set E [expr $SE + $EE]
   set lastEE $EE
   set lastSE $SE
   set E0 $E
   puts "CFAFLEXMC) E0 $E0"
   set keep_cycling 1
   if { $EE < $sstop && $SE < $dstop } {
      set keep_cycling 0
   }
   set nrb [bondstruct_getnrb $bs]
   #set tnacc 0
   #set ngc 0
   set nacc 0
   for {set cyc 0} { $cyc < $maxcycles && $keep_cycling == 1 } { incr cyc } {
      set SAVEPOS [$rotsel get {x y z}] ; # modify so only atoms on right-side of this bond are saved
      puts -nonewline "rotating..."
      flush stdout
      set start [clock microseconds]
      for {set r 0} {$r < $nrb} {incr r} {
         #get a random active bond
        #set rb [irand_dom 0 [expr $nrb-1]]
        # get a random rotation angle
        set av [expr $maxanglestep * [irand_dom -$anglegradationfactor $anglegradationfactor]]
        # get this active bonds index in the bondstruct (why don't I just delete nonrotatable bonds?)
#        set rr [bondstruct_r2b $bs $rb]
        set rr [bondstruct_r2b $bs $r]
        # execute the rotation
        #if { [bondstruct_isactive $bs $rr] } {
           bondrot_by_index $bs $molid $rr $av
        #}
      }
      set stop [clock microseconds]
      puts "...done [expr ($stop-$start)/1.e6] s"
      set SE 0.0
      foreach i $ilist j $jlist {
        if { $i != $j } {
          set dattr [measure bond [list $i $j]]
          set SE [expr $SE+0.5*$k*pow($dattr,2)]
        }
      }
      puts -nonewline "roughenergy..."
      flush stdout
      set start [clock microseconds]
      set EE [roughenergy $rotnoh $ljcutoff $ljsigma $ljepsilon $ljshift $bs $ls]
      set stop [clock microseconds]
      puts "...done [expr ($stop-$start)/1.e6] s"
      set E [expr $SE + $EE]
      # do metropolis
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
        $rotsel set {x y z} $SAVEPOS
      } else {
          # accept the move
          set E0 $E
          puts "CFAFLEXMC) cyc $cyc na $nacc [format "ar=%.5f" [expr (1.0*$nacc)/($cyc+1)]] [format "attr-pnlty= %.2f " $SE] [format "strc-pnlty=%.2f" $EE]"
          if { [expr $nacc % $logevery == 0 ] && $logid != -1 } {
            log_addframe $molid $logid
            if { [expr $nacc % $logsaveevery == 0] } {
                set loga [atomselect $logid all]
                animate write dcd "tmp.dcd" waitfor all sel $loga $logid
                $loga delete
            }
          }
      }
      set lastEE $EE
      set lastSE $SE
      if { $EE < $sstop && $SE < $dstop } {
        set keep_cycling 0
      }
   }
   puts -nonewline "CFAFLEXMC) Stop at cycle $cyc: "
   if { $i != $j } {
     puts -nonewline "attr dst: [format "%.2f" [measure bond [list $i $j]]] [format "attr-pnlty %.2f " $lastSE]"
   }
   puts "[format "strc-pnlty %.2f" $lastEE]"
   free_bondstruct $bs
   roughenergy_cleanup $ls
}

proc ladd {l} {
    set total 0.0
    foreach nxt $l {
        set total [expr $total + $nxt]
    }
    return $total
}

proc check_pierced_rings { molid ringsize TOL } {
  set r6 [atomselect $molid "ringsize $ringsize from all"]
  set r6i [$r6 get index]
  set i 0
  foreach ii $r6i {
    set r6o($ii) $i
    incr i
  }
  set r6x [$r6 get x]
  set r6y [$r6 get y]
  set r6z [$r6 get z]

  for { set ri 0 } { $ri < [llength $r6i] } { incr ri $ringsize } {
    #puts "ring $ri"
    set this_ri {}
    set this_rx {}
    set this_ry {}
    set this_rz {}
    for { set t 0 } { $t < $ringsize } { incr t } {
      lappend this_ri [lindex $r6i [expr $ri + $t]]
      lappend this_rx [lindex $r6x [expr $ri + $t]]
      lappend this_ry [lindex $r6y [expr $ri + $t]]
      lappend this_rz [lindex $r6z [expr $ri + $t]]
    }
    #puts "this_ri $this_ri"
    set this_rr {}
    foreach x $this_rx y $this_ry z $this_rz {
      lappend this_rr [list $x $y $z]
    }
    #puts "this_rr $this_rr"
    set this_com [list [ladd $this_rx] [ladd $this_ry] [ladd $this_rz]]
    set this_com [vecscale $this_com [expr 1./$ringsize]]
    set this_b12 [vecsub [lindex $this_rr 0] [lindex $this_rr 1]]
    set this_b23 [vecsub [lindex $this_rr 1] [lindex $this_rr 3]]
    set c123 [veccross $this_b12 $this_b23]
    set lc123 [veclength $c123]
    set chat123 [vecscale $c123 [expr 1.0/$lc123]]
    #puts "ring $this_ri : $this_com : $chat123"
    set neigh [atomselect $molid "(protein or glycan) and same residue as (exwithin 4.0 of index $this_ri)"]
    set nb [$neigh getbonds]
    set na [$neigh get index]
    set nax [$neigh get x]
    set nay [$neigh get y]
    set naz [$neigh get z]
    set i 0
    foreach a $na {
      set ord($a) $i
      incr i
    }
    set ln 0
    set ndots 0
    # fix to exclude atoms in the ring from bond definition!
    foreach a $na bl $nb {
      if { [lsearch $this_ri $a] !=  -1 } {
        continue
      }
      if { [expr $ln%100 == 0] } {
        puts -nonewline "."
        incr ndots
        if { $ndots > 80 } {
          puts ""
          set ndots 0
        }
        flush stdout
      }
      incr ln
      set ai $ord($a)
      set apos [list [lindex $nax $ai] [lindex $nay $ai] [lindex $naz $ai]]
      foreach b $bl {
        if { [lsearch $this_ri $b] != -1 } {
          continue
        }
        if { [ expr $b < $a ] } {
          continue
        }
        if { [lsearch $na $b] != -1 } {
          set bi $ord($b)
          set bpos [list [lindex $nax $bi] [lindex $nay $bi] [lindex $naz $bi]]
          set avpos [vecscale [vecadd $apos $bpos] 0.5]
          set avec [vecsub $this_com $apos]
          set bvec [vecsub $this_com $bpos]
          set crit1 [expr [veclength [vecsub $avpos $this_com]] < $TOL]
          set d1 [vecdot $avec $chat123]
          set d2 [vecdot $bvec $chat123]
          set crit2 [expr ($d1*$d2)<0]
          if { $crit1 && $crit2 } {
            puts ""
            set piercer [atomselect $molid "index $a $b"]
            set b_resnames [$piercer get resname]
            set b_names [$piercer get name]
            set b_segnames [$piercer get segname]
            set b_resids [$piercer get resid]
            set b_chain [$piercer get chain]
            set piercee [atomselect $molid "index $this_ri"]
            set first "[lindex $b_chain 0]([lindex $b_segnames 0])_[lindex $b_resnames 0][lindex $b_resids 0][lindex $b_names 0]($a)"
            set second "[lindex $b_chain 1]([lindex $b_segnames 1])_[lindex $b_resnames 1][lindex $b_resids 1][lindex $b_names 1]($b)"
            set ring "[lindex [$piercee get chain] 0]([lindex [$piercee get segname] 0])_[lindex [$piercee get resname] 0][lindex [$piercee get resid] 0]"
            puts "Bond $first-$second pierces $ring"
          }
        }
      }
    }
    $neigh delete
  }
}

proc difference { a b } {
  set r {}
  foreach i $a {
    if { [lsearch $b $i] == -1 } {
      lappend r $i
    }
  }
  return $r
}

proc ligateCN { molid residueC residueN } {
  set jsel [atomselect $molid "(residue $residueC and name C OT1 OT2) or (residue $residueN and name N HT1 HT2 HT3)"]
  set an [$jsel get name]
  set segnames [$jsel get segname]
  set resids [$jsel get resid]
  for { set i 0 } { $i < [llength $an] } { incr i } {
    set index([lindex $an $i]) $i
    set resid([lindex $an $i]) [lindex $resids $i]
    set segname([lindex $an $i]) [lindex $segnames $i]
  }

  # pick the one OT and the one HT that would give the most "trans" peptide bond
  set dimax 0.0
  set thetwo {}
  foreach o { OT1 OT2 } {
    foreach h { HT1 HT2 HT3 } {
      set thisdi [expr abs([measure dihed [list $index($o) $index(C) $index(N) $index($h)]])]
      if { $thisdi > $dimax } {
        set dimax $thisdi
        set thetwo [list $o $h]
      }
    }
  }
  set fullset { OT1 OT2 HT1 HT2 HT3 }
  # symmetric difference
  set deleteus [difference $fullset $thetwo]
  puts "thetwo $thetwo deleteus $deleteus"
  foreach db $deleteus {
    delatom $segname($db) $resid($db) $db
  }

}