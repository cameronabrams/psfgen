# VMD/psfgen script for generating psf/pdb pair for a single BMS-529
# molecule extracted from the 5u7o PDB entry
#
# cameron f abrams (c) 2018
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

mol new 5u7o.pdb
set a [atomselect top "resname 83J"]
#$a set resname "SUCR"
$a set chain X
$a set resid 1
$a set resname B529 
$a moveby [vecscale -1 [measure center $a weight mass]]

$a writepdb "b529.pdb"
lappend LOCALFILES b529.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology $PSFGEN_BASEDIR/charmm/bms529.str

segment X {
   pdb b529.pdb
}
coordpdb b529.pdb X
#guesscoord
writepsf "my_b529.psf"
writepdb "my_b529.pdb"

foreach f $LOCALFILES {
  exec /bin/rm -f $f
}


exit

