### globglob.tcl
### psfgen script for generating a single fusion of two globular proteins
###
### Cameron F. Abrams cfa22@drexel.edu

# check for file existence and die if not found
proc filechk { fn } {
    if { ! [file exists $fn ] } {
        puts "ERROR: $fn not found."
        exit 1
    }
}

# give me a random number between min and max
proc irand_dom { min max } { 
   return [expr int(rand() * ($max - $min + 1)) + $min]
}

# rotate i->j bond of $sel of $molid by $deg (only atoms j and beyond move)
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

# end-to-end distance of sel (alpha carbon reference)
proc e2e { molid sel } {
    set resids [lsort -unique -integer [$sel get resid]]
    set nt [atomselect $molid "[$sel text] and resid [lindex $resids 0] and name CA"]  
    set ct [atomselect $molid "[$sel text] and resid [lindex $resids end] and name CA"]
    set np [lindex [$nt get {x y z}] 0]
    set cp [lindex [$ct get {x y z}] 0]
    return [veclength [vecsub $np $cp]]
}

# randomly assign phi-psi angles of all residues in $sel; no overlaps; 
# end-to-end must be greater than mine2e
proc random_loop { molid sel mine2e } {
   set ni {}
   set ci {}
   set cai {}
   set nmaxtry 100
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
#     puts "$psi $phi"
     flush stdout
     set PSIM [measure dihed $psi]
     set PHIM [measure dihed $phi]
#     puts "psi $psi ([format %.2f $PSIM]) phi $phi ([format %.2f $PHIM])"
     # pick a random (phi,psi) from ramachandran plot
     set clashing 1
     set ntry 0
     while { $clashing == 1 && $ntry < $nmaxtry } {
         set clashing 0
         set PHIR [irand_dom -340 -60]
         set PSIR [irand_dom -60 300]
         set SPHI [expr ($PHIR) - ($PHIM)]
         set SPSI [expr ($PSIR) - ($PSIM)]
         genbondrot $molid $sel [lindex $psi 1] [lindex $psi 2] $SPSI
         genbondrot $molid $sel [lindex $phi 1] [lindex $phi 2] $SPHI
         # if clashes are produced, try again
         set contacts [measure contacts 1.5 $sel]
         set coni [lindex $contacts 0]
         set conj [lindex $contacts 1]
         if { [llength $coni] > 0 } {
             set clashing 1
             puts "Clash detected...trying again..."
         }
         if { [e2e $molid $sel] < $mine2e } {
             set clashing 1
             puts "Ends too close...trying again..."
         }
         incr ntry
     }
     if { $ntry == $nmaxtry } {
         puts "ERROR: could not produce a viable configuration after $nmaxtry tries"
         exit 1
     } else {
         puts "Set angles of resid $in"
     }

#     puts "$SPHI $SPSI"
   }
}

package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

# first arg is pdb for first globule
set glob1pdb [lindex $argv 0]
# second arg is pdb for second globule
set glob2pdb [lindex $argv 1]

filechk $glob1pdb
filechk $glob2pdb

mol new $glob1pdb
set glob1id 0
set g1 [atomselect $glob1id "protein and chain A"]
$g1 num
# save chain A of glob1
$g1 writepdb "glob1A.pdb"
# make list of resids for glob1
set resids1 [lsort -unique -integer [$g1 get resid]]

mol new $glob2pdb
set glob2id 1
set g2 [atomselect $glob2id "protein and chain A"]
$g2 num
# save chain A of glob2
$g2 writepdb "glob2A.pdb"
# make list of resids for glob2
set resids2 [lsort -unique -integer [$g2 get resid]]

# add 51 alanines to the C-terminus of glob1
segment A {
    pdb glob1A.pdb
    for { set i 0 } { $i < 51 } { incr i } {
        residue [expr $i + [lindex $resids1 end] + 1] ALA A
    }
}
coordpdb glob1A.pdb A

guesscoord
writepsf "try1.psf"
writepdb "try1.pdb"

resetpsf
mol delete $glob1id
mol new try1.pdb
set glob1id [molinfo top get id]
# select the linker
set linker [atomselect $glob1id "resid [expr [lindex $resids1 end] + 1] to [expr [lindex $resids1 end] + 52]"]
# fiddle with all the phi-psi angles
random_loop $glob1id $linker 20.0
# save it
[atomselect $glob1id "all"] writepdb "try2.pdb"

mol delete $glob1id
resetpsf
mol new "try2.pdb"
set glob1id [molinfo top get id]
set g1 [atomselect $glob1id "protein and chain A"]
set resids1 [lsort -unique -integer [$g1 get resid]]

# shift resids of glob2 residues
set newg2n [expr [lindex $resids1 end]+1]
set newresids2 {}
set allresids2 [$g2 get resid]
foreach r $allresids2 {
    lappend newresids2 [expr $r + $newg2n]
}
$g2 set resid $newresids2

# shift position of 
# locate the N of first res of so that it overlaps one of the two OT's on cterm
set nnter2 [atomselect $glob2id "resid [lindex $newresids2 0] and name N"]
puts "nnter2 [$nnter2 num]"
set ccter1 [atomselect $glob1id "resid [lindex $resids1 end] and name OT1"]
puts "ccter1 [$ccter1 num]"
set nr [lindex [$nnter2 get {x y z}] 0]
set or [lindex [$ccter1 get {x y z}] 0]
puts "$nr $or"
set shift [vecsub $or $nr]
puts "Moving glob2 by $shift..."
$g2 moveby $shift
set c2 [$g2 measure center]
$g2 writepdb "try3.pdb"
resetpsf

segment A {
    pdb try2.pdb
    pdb try3.pdb
}
coordpdb try2.pdb A
coordpdb try3.pdb A
guesscoord
writepsf "fusion.psf"
writepdb "fusion.pdb"

exit

