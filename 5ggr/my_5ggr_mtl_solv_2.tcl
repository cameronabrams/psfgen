package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd.str

readpsf my_5ggr.psf
coordpdb prot.pdb
pdbalias residue DMAN DMANOL
mol new my_5ggr_mtl.pdb ; # output of packmol
set segs { M WA WB WC WD WE WF WG WI WJ WK WN  I }
set chns { M  A  B  C  D  E  F  G  I  J  K  N  O }
foreach sn $segs cn $chns {
  set s [atomselect top "segname $sn"] 
  $s set chain $cn
  $s writepdb "${sn}.pdb"
  segment $sn {
    if { [string index $sn 0] == "W" } {
      auto none
    }
    pdb ${sn}.pdb
  }
  coordpdb ${sn}.pdb $sn
}

regenerate angles dihedrals

set outputname my_5ggr_mtl_solv
writepsf "${outputname}.psf"
writepdb "${outputname}.pdb"

# generate a special PDB with certain residue alpha-carbons tagged for use
# by a colvars module input file for orientational and positional restraints
mol new ${outputname}.psf
mol addfile ${outputname}.pdb
set a [atomselect top all]
$a set beta 0.0
set b [atomselect top "name CA and ((chain H and resid 2 to 128 134 to 213) or (chain L and resid 1 to 212))"]
$b set beta 1.0
$b writepdb ${outputname}_caB1.pdb

quit
