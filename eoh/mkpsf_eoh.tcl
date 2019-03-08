# VMD/psfgen script for generating psf/pdb pair for a single ethanol
# molecule extracted from the 3tod PDB entry
#
# cameron f abrams (c) 2019
# cfa22@drexel.edu
# drexel university
# chemical and biological engineering
#
# check for base directory name variable;
# if not set, use default
if {![info exists PSFGEN_BASEDIR]} {
  # see if user set an environment variable
  if {[info exists env(PSFGEN_BASEDIR)]} {
      set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
  } else {
      set PSFGEN_BASEDIR $env(HOME)/research/psfgen
  }
}
set LOCALFILES {}

mol new 3tod.pdb
set a [atomselect top "resname EOH and resid 701"]
#$a set resname "SUCR"
$a set chain E
$a set resid 1
$a set resname ETOH 
$a moveby [vecscale -1 [measure center $a weight mass]]

$a writepdb "eoh.pdb"
lappend LOCALFILES eoh.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_lipid.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf

segment E {
   pdb eoh.pdb
}
coordpdb eoh.pdb E
guesscoord
writepsf "my_eoh.psf"
writepdb "my_eoh.pdb"

foreach f $LOCALFILES {
  exec /bin/rm -f $f
}


exit

