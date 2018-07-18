# VMD/psfgen script for generating psf/pdb pair for a single 
# BNM-III-170 molecule extracted from the 5f4p PDB entry
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

mol new 5f4p.pdb
set a [atomselect top "resname 5VG"]
$a set chain X
$a set resid 1
$a set resname BNM3
$a moveby [vecscale -1 [measure center $a weight mass]]

$a writepdb "bnm.pdb"
lappend LOCALFILES bnm.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology ${PSFGEN_BASEDIR}/charmm/bnm.str

foreach nbad [exec grep 5VG 5f4p.pdb | grep -w "5VG A" | grep HETATM | cut -b 13-16] ngood [exec grep ATOM ${PSFGEN_BASEDIR}/charmm/bnm.str | cut -b 6-9 | grep -v ^H] {
  pdbalias atom BNM3 $nbad $ngood
} 

segment X {
   pdb bnm.pdb
}
coordpdb bnm.pdb X
guesscoord
writepsf "my_bnm.psf"
writepdb "my_bnm.pdb"

foreach f $LOCALFILES {
  exec /bin/rm -f $f
}


exit

