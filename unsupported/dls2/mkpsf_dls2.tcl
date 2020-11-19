# VMD/psfgen script for generating psf/pdb pair for a single DLS2
# molecule -- repeat unit for linker polymer in BMS-linker-Trp3 DAVEI
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

mol new dls2_raw.pdb
set a [atomselect top "all"]
$a set chain X
$a set resid 3
$a set resname DLS2 
$a moveby [vecscale -1 [measure center $a weight mass]]

$a writepdb "dls2.pdb"
#lappend LOCALFILES dls2.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology $PSFGEN_BASEDIR/charmm/dls2.str

segment X {
   pdb dls2.pdb
}
coordpdb dls2.pdb X
guesscoord

writepsf "my_dls2.psf"
writepdb "my_dls2.pdb"

foreach f $LOCALFILES {
  exec /bin/rm -f $f
}


exit

