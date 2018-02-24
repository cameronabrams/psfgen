# VMD/psfgen script for generating psf/pdb pair for a single D-mannitol
# molecule extracted from the 1m2w PDB entry
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

mol new 1m2w.pdb
set a [atomselect top "resname MTL and chain A"]
#$a set resname "SUCR"
$a set chain M
$a set resid 1
$a set resname DMTL 
$a moveby [vecscale -1 [measure center $a weight mass]]

$a writepdb "mtl.pdb"
lappend LOCALFILES mtl.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
pdbalias residue DMTL DMANOL

segment M {
   pdb mtl.pdb
}
coordpdb mtl.pdb M
guesscoord
writepsf "my_mtl.psf"
writepdb "my_mtl.pdb"

foreach f $LOCALFILES {
  exec /bin/rm -f $f
}


exit

