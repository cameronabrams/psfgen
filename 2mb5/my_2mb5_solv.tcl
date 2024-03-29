# VMD script for generating a solvated system given the protein psf/pdb
#
# cameron f abrams (c) 2017
# drexel university chemical and biological engineering
#
#
set pad 10 ; #  pad in angstroms

set inputname my_2mb5
set PSF ${inputname}.psf
set outputname1 ${inputname}_wb
set outputname2 ${inputname}_i

mol new $PSF
mol addfile ${inputname}_vac.coor

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

$a writepdb ${inputname}.pdb

package require solvate
package require autoionize

solvate $PSF ${inputname}.pdb -minmax $box -o ${outputname1}
autoionize -psf ${outputname1}.psf -pdb ${outputname1}.pdb -neutralize -o ${outputname2}

# generate an input file for the first solvated MD simulation
# namd config file
set fp [open "cell.inp" "w"]
puts $fp "cellbasisvector1 [lindex $basisvec 0] 0 0"
puts $fp "cellbasisvector2 0 [lindex $basisvec 1] 0"
puts $fp "cellbasisvector3 0 0 [lindex $basisvec 2]"
puts $fp "cellorigin $origin"
close $fp
puts "Generated cell.inp."

# generate a special PDB with certain reside alpha-carbons tagged for use
# by a colvars module input file for orientational and positional restraints
mol new ${outputname2}.psf
mol addfile ${outputname2}.pdb
set a [atomselect top all]
$a set beta 0.0
set b [atomselect top "name CA and chain A and resid 1 to 153"] 
$b set beta 1.0
$b writepdb ${inputname}_caB1.pdb
quit

