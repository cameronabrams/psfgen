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

# computes overlap energy between atoms in sel1 and sel2.  The "my_roughenergy" function (implemented in C)
# uses a repulsive pair potential of the form A*(x-cut)^2 for x<cut.  Residue index lists are
# sent so the my_roughenergy does not compute pair interactions for atoms in the same residue 
proc roughenergy { sel1 sel2 cut }  {
  set E 0.0
  if { [$sel1 num] > 0 && [$sel2 num] > 0 } {
   set r1 [intListToArray [$sel1 get residue]]
   set r2 [intListToArray [$sel2 get residue]]
   set x1 [ListToArray [$sel1 get x]]
   set x2 [ListToArray [$sel2 get x]]
   set y1 [ListToArray [$sel1 get y]]
   set y2 [ListToArray [$sel2 get y]]
   set z1 [ListToArray [$sel1 get z]]
   set z2 [ListToArray [$sel2 get z]]
   set E [my_roughenergy $r1 $x1 $y1 $z1 [$sel1 num] $r2 $x2 $y2 $z2 [$sel2 num] $cut]
   delete_arrayint $r1
   delete_arrayint $r2
   delete_array $x1
   delete_array $x2
   delete_array $y1
   delete_array $y2
   delete_array $z1
   delete_array $z2
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
   set p2 [lindex [$cg get {x y z}] 0]
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
     [atomselect $molid all] writepdb "tmp.pdb"
     animate read pdb tmp.pdb $logid
     puts "Molid $molid - logging molecule $logid has [molinfo $logid get numframes] frames."
     exec rm -f tmp.pdb
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
proc do_loop_mc { residueList c molid k r0 env rcut maxcycles temperature iseed logid } {

  set msel [atomselect $molid "chain $c and residue $residueList"]
  set mselnoh [atomselect $molid "chain $c and residue $residueList and noh"]
  set mselca [atomselect $molid "chain $c and residue $residueList and name CA"]
  set envca [atomselect $molid "([$env text]) and name CA"]

  set rend [lindex $residueList end]
  set nres [llength $residueList]

  set co [atomselect $molid "chain $c and residue $rend and name C"]
  set ca [atomselect $molid "chain $c and residue $rend and name CA"]
  set idx [list [$ca get index] [$co get index]]
  set ibl [format "%.2f" [measure bond $idx]]
  puts "CFALOOPMC) residue $rend CA_[$ca get index]-C_[$co get index]: $ibl A"

  expr srand($iseed)

  set nacc 0

  set SE [expr 0.5*$k*pow([measure bond $idx]-$r0,2)]
  #set EE [roughenergy $mselnoh $env $rcut]
  set EE [roughenergy $mselca $envca $rcut]
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
    #set EE [roughenergy $mselnoh $env $rcut]
    set EE [roughenergy $mselca $envca $rcut]
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
      log_addframe $molid $logid
    }
  }
  puts "CFALOOPMC) ($rend) cyc $cyc na $nacc ([format "%.5f" [expr (1.0*$nacc)/($cyc+1)]]) CA-C: [format "%.2f" [measure bond $idx]]"
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
proc do_flex_mc { molid msel ri rj fa k i j envsel rcut maxcycles temperature iseed logid } {

   set bl [$msel getbonds]
   set il [$msel get index]
   set bs [make_bondstruct $molid $msel $ri $rj]
   bondstruct_deactivate_by_fixed $bs $fa
#   print_bondlist $bs

#   puts "ri $ri"
#   puts "rj $rj"
#   puts "il [llength $il] : $il"
#   puts "bl [llength $bl] : $bl"
   
   if { $i != $j } { 
     puts "CFAFLEXMC) Initial attractor distance [format "%.2f" [measure bond [list $i $j]]] A"
   }

   flush stdout

   expr srand($iseed)

   set nacc 0

   set SE 0.0
   if { $i != $j } {
     set SE [expr 0.5*$k*pow([measure bond [list $i $j]],2)]
   }
   set EE [roughenergy $msel $envsel $rcut]
   set E [expr $SE + $EE]
   set E0 $E
#   puts "CFAFLEXMC) E0 $E0"
   for {set cyc 0} { $cyc < $maxcycles } { incr cyc } {
      # save coordinates
      set SAVEPOS [$msel get {x y z}]
      set nrot 0
      for {set r 0} {$r < [bondstruct_getnb $bs] } {incr r} {
         #set av [expr 60 * [irand_dom 1 5]]
         set av [expr 6 * [irand_dom -5 5]]
        # puts "cyc $cyc bond $r deg $av"
         if { [bondstruct_isactive $bs $r] } {
           bondrot_by_index $bs $molid $r $av
           set nrot [expr $nrot + 1]
         }
      }
    
      if { $nrot == 0 } {
         puts "ERROR: no rotations performed"
         exit
      }
      if { $i != $j } {
        set SE [expr 0.5*$k*pow([measure bond [list $i $j]],2)]
      } else {
        set SE 0.0
      }
      set EE [roughenergy $msel $envsel $rcut]
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
        log_addframe $molid $logid
      }
   }
   if { $i != $j } {
     puts -nonewline "attr dst: [format "%.2f" [measure bond [list $i $j]]] [format "attr-pnlty %.2f " $SE]"
   }
   puts "[format "strc-pnlty %.2f" $EE]"
   free_bondstruct $bs
}

#proc check_pierced_rings molid {
  # molid is a molecule assumed to have some residues with rings and perhaps glycans

  # this will search the list of bonds and for each, it will search all rings within 5.0 
  # angstroms of the bond to determine if the bond pierces one of those rings.

#  set a [atomselect $molid "protein or glycan"]
#  set ai [$a get index]
#  set bl [$a get bonds]
# foreach i $ai b $bl {
#    set r5 [atomselect $molid "(same residue as within 5 of index $i) and ringsize 5"]
#    if {[expr [$r5 num] > 0]} {
# MUCH WORK IS TO DO HERE
#    }
#  }
#}

proc glycan_rotatables { molid selstr } {
  set a [atomselect $molid "$selstr"]
  set bl [$a getbonds]
  set at [$a get index]
  set r6 [atomselect $molid "ringsize 6 from ($selstr)"]
  set ri6 [$r6 get index]
  puts "ringsize sel has [$r6 num] atoms"
  set lefts [list]
  set rights [list]
  foreach i $at b $bl {
    foreach j $b {
      if { [expr $j > $i] } {
        puts "considering $i - $j"
        set inring 0
        for { set r 0 } { $r < [llength $ri6] } { incr r 6 } {
            #puts "   considering ring [expr $r/6]"
            set i_in 0
            set j_in 0
            for { set ri $r } { $ri < [expr $r + 6] } { incr ri } {
               if { $i == [lindex $ri6 $ri] } {
                  set i_in 1
               }
               if { $j == [lindex $ri6 $ri] } {
                  set j_in 1
               }
            }
            if { $i_in == 1 && $j_in == 1 } {
              #puts "bond $i - $j lives in ring [expr $r/6]"
              set inring 1
              break
            }
        }
        if { $inring == 1 } {
           # puts "$i $j is a ring edge"
        } elseif { $i_in == 0 } {
            #puts "$j is a ring vertex but $i is not"
            lappend lefts $j
            lappend rights $i
        } elseif { $j_in == 0 } {
            #puts "$i is a ring vertex but $j is not"
            lappend lefts $i
            lappend rights $j
        }
      }
    }
  }
  return [list $lefts $rights]
}
