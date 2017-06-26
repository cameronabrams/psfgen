# VMD script for generating a solvated system given the protein psf/pdb
#
# cameron f abrams (c) 2017
# drexel university chemical and biological engineering
#
#
set pad 10 ; #  pad in angstroms

mol new my_3tgq.psf
mol addfile my_3tgq_relax.coor

set a [atomselect top all]

set c [measure center $a]

$a moveby [vecscale $c -1]

set box { { ? ? ? } { ? ? ? } }
set basisvec { ? ? ? }
set origin { ? ? ? }
set minmax [measure minmax $a]

foreach d {0 1 2} {
  lset box 0 $d [expr [lindex $minmax 0 $d] - $pad]
  lset box 1 $d [expr [lindex $minmax 1 $d] + $pad]
  lset basisvec $d [expr [lindex $box 1 $d ] - [lindex $box 0 $d]] 
  lset origin $d [expr 0.5*([lindex $box 1 $d ] + [lindex $box 0 $d])] 
}

$a writepdb my_3tgq_relax.pdb

package require solvate
package require autoionize

solvate my_3tgq.psf my_3tgq_relax.pdb -minmax $box -o my_3tgq_wb
autoionize -psf my_3tgq_wb.psf -pdb my_3tgq_wb.pdb -neutralize -o my_3tgq_i

# generate an input file for the first solvated MD simulation
# namd config file
set fp [open "cell.inp" "w"]
puts $fp "cellbasisvector1 [lindex $basisvec 0] 0 0"
puts $fp "cellbasisvector2 0 [lindex $basisvec 1] 0"
puts $fp "cellbasisvector3 0 0 [lindex $basisvec 2]"
puts $fp "cellorigin $origin"
close $fp
puts "Generated cell.inp."

quit

