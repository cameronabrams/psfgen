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

mol new dls1_raw.pdb
set a [atomselect top "all"]
#$a set resname "SUCR"
$a set chain X
$a set resid 2
$a set resname DLS1 
$a moveby [vecscale -1 [measure center $a weight mass]]

$a writepdb "dls1.pdb"
lappend LOCALFILES dls1.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology $PSFGEN_BASEDIR/charmm/dls1.str

segment X {
   pdb dls1.pdb
}
coordpdb dls1.pdb X
guesscoord
writepsf "my_dls1.psf"
writepdb "my_dls1.pdb"

foreach f $LOCALFILES {
  exec /bin/rm -f $f
}


exit

