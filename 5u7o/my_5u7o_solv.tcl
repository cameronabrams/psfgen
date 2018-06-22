# VMD script for generating solvated/neutralized system from the protein psf/pdb
#
# cameron f abrams (c) 2017-2018
# cfa22@drexel.edu
# drexel university
# chemical and biological engineering

set pad 10; # pad in angstroms

set inputname my_5u7o
set PSF ${inputname}.psf
set outputname1 ${inputname}_wb
set outputname2 ${inputname}_i

mol new $PSF
mol addfile ${inputname}_vac.coor

set a [atomselect top all]

set c [measure center $a]

$a moveby [vecscale $c -1.0]

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

# generate a special PDB with certain residue alpha-carbons tagged for use
# by a colvars module input file for orientational and positional restraints
mol new ${outputname2}.psf
mol addfile ${outputname2}.pdb
set a [atomselect top all]
$a set beta 0.0
set b [atomselect top "name CA and chain G K R and resid 253 to 301 322 to 396 411 to 475"] 
set rcm [measure center $b]
$b set beta 1.0
$b writepdb "my_5u7o_caB1.pdb"

set fp [open "my_5u7o_colvars_op.inp" "w"]
puts $fp "colvarstrajfrequency 100
colvar {
  name prot_orientation
  orientation {
    atoms {
      psfSegID G G G K K K R R R
      atomNameResidueRange CA 253-301
      atomNameResidueRange CA 322-396
      atomNameResidueRange CA 411-475
      atomNameResidueRange CA 253-301
      atomNameResidueRange CA 322-396
      atomNameResidueRange CA 411-475
      atomNameResidueRange CA 253-301
      atomNameResidueRange CA 322-396
      atomNameResidueRange CA 411-475
    }
    refpositionsfile my_5u7o_caB1.pdb
    refpositionscol B
  }
}

colvar {
  name prot_position
    distance {
        group1 {
           psfSegID G G G K K K R R R
           atomNameResidueRange CA 253-301
           atomNameResidueRange CA 322-396
           atomNameResidueRange CA 411-475
           atomNameResidueRange CA 253-301
           atomNameResidueRange CA 322-396
           atomNameResidueRange CA 411-475
           atomNameResidueRange CA 253-301
           atomNameResidueRange CA 322-396
           atomNameResidueRange CA 411-475
        }
        group2 {
            dummyatom (0.0, 0.0, [lindex $rcm 2])
        }
    }
}
harmonic {
  name prot_fixor
  colvars prot_orientation
  forceconstant 100.0
  centers (1.0 , 0.0 , 0.0 , 0.0)
}

harmonic {
    name prot_nodrift
    colvars prot_position
    forceconstant 100.0
    centers 0.0
}"
close $fp

quit
