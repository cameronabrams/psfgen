# VMD script for generating solvated/neutralized system from the protein psf/pdb
#
# cameron f abrams (c) 2017-2020
# cfa22@drexel.edu
# drexel university
# chemical and biological engineering

set scriptname $argv0

set pad 10; # pad in angstroms
set pdb "empty.pdb"
set psf "empty.psf"
set outpre "ionized"
set cubic 0

for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pad" } {
       incr i
       set pad [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-cubic" } {
       set cubic 1
    }
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set pdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-psf"} {
       incr i
       set psf [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-outpre"} {
       incr i
       set outpre [lindex $argv $i]
    }
}

set outputname1 ${outpre}_wb

mol new $psf
mol addfile $pdb

set a [atomselect top all]

set box { { ? ? ? } { ? ? ? } }
set basisvec { ? ? ? }
set origin { ? ? ? }

set minmax [measure minmax $a]
set maxspan [expr -9999]
foreach d {0 1 2} {
   set thisspan [expr [lindex $minmax 1 $d ] - [lindex $minmax 0 $d]]
   if { $thisspan > $maxspan } {
      set maxspan $thisspan
   }
}
vmdcon -info "${scriptname}: Maximum span $maxspan Angstroms"

set sympad [list 0 0 0]
if { $cubic == 1 } {
   vmdcon -info "${scriptname}: Enforcing cubic box"
   foreach d {0 1 2} {
      set thisspan [expr [lindex $minmax 1 $d ] - [lindex $minmax 0 $d]]
      lset sympad $d [expr 0.5*($maxspan-$thisspan)]
   }
}
foreach d {0 1 2} {
  lset box 0 $d [expr [lindex $minmax 0 $d] - $pad - [lindex $sympad $d]]
  lset box 1 $d [expr [lindex $minmax 1 $d] + $pad + [lindex $sympad $d]]
  lset basisvec $d [expr [lindex $box 1 $d ] - [lindex $box 0 $d]] 
  lset origin $d [expr 0.5*([lindex $box 1 $d ] + [lindex $box 0 $d])] 
}

package require solvate
package require autoionize
psfcontext mixedcase

solvate $psf $pdb -minmax $box -o ${outputname1}
autoionize -psf ${outputname1}.psf -pdb ${outputname1}.pdb -neutralize -o ${outpre}

# generate an input file for the first solvated MD simulation
# namd config file
set fp [open "cell.inp" "w"]
puts $fp "cellbasisvector1 [lindex $basisvec 0] 0 0"
puts $fp "cellbasisvector2 0 [lindex $basisvec 1] 0"
puts $fp "cellbasisvector3 0 0 [lindex $basisvec 2]"
puts $fp "cellorigin $origin"
close $fp
vmdcon -info "${scriptname}: Generated cell.inp."

quit
