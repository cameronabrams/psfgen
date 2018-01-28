# VMD/psfgen script for generating psf/pdb pair for a single sucrose
# molecule extracted from the 5o8l PDB entry
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

mol new 5o8l.pdb
set a [atomselect top "resname SUC"]
#$a set resname "SUCR"
$a set chain S
$a set resid 1
$a moveby [vecscale -1 [measure center $a weight mass]]

set aglc [atomselect top "name O1 C1 C2 O2 C3 O3 C4 O4 C5 O5 C6 O6"]
$aglc set resid 1
$aglc set resname AGLC
set bfru [atomselect top "name C1' O1' C2' C3' O3' C4' O4' C5' C6' O6' O2'"]
set name [$bfru get name]
set newname {}
foreach n $name {
   lappend newname [string trim $n "'"]
}

$bfru set name $newname
$bfru set resid 2
$bfru set resname BFRU
[atomselect top "resname BFRU and name O2"] set name "O5"

$aglc writepdb "aglc.pdb"
$bfru writepdb "bfru.pdb"
lappend LOCALFILES aglc.pdb
lappend LOCALFILES bfru.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_carb.rtf

segment A {
   pdb aglc.pdb
   pdb bfru.pdb
}
coordpdb aglc.pdb A
coordpdb bfru.pdb A

patch SUCR A:1 A:2
guesscoord
writepsf "my_sucr.psf"
writepdb "my_sucr.pdb"

foreach f $LOCALFILES {
  exec /bin/rm -f $f
}


exit

