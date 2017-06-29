# VMD/psfgen script for generating psf/pdb pair for PDB 1f7a
# dimeric HIV-1 protease with a bound substrate mimic
# and a deactivating mutation reverted to wild-type
#
# cameron f abrams (c) 2017
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
# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/tcl/loopmc.tcl
set LOCALFILES {}

mol new 1f7a.pdb

set all [atomselect top "all"]
set A [atomselect top "protein and chain A"]
$A writepdb A.pdb; lappend LOCALFILES A.pdb
set B [atomselect top "protein and chain B"]
$B writepdb B.pdb; lappend LOCALFILES B.pdb
set P [atomselect top "protein and chain P"]
$P writepdb P.pdb; lappend LOCALFILES P.pdb
set W [atomselect top "water and within 3.0 of protein"]
$W set name OH2
$W set resname TIP3
$W writepdb W.pdb; lappend LOCALFILES W.pdb

package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $PSFGEN_BASEDIR/charmm/prod.top
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str

pdbalias atom ILE CD1 CD
pdbalias residue HIS HSD

segment A { 
  pdb A.pdb
  mutate 25 ASP
}
coordpdb A.pdb A

segment B {
  pdb B.pdb
  mutate 25 ASP
}
coordpdb B.pdb B

segment P {
  pdb P.pdb
}
coordpdb P.pdb P

segment W {
  pdb W.pdb
}
coordpdb W.pdb W

patch ASPP A:25
patch ASPP B:25
patch PROD A:44
patch PROD B:44

guesscoord

regenerate angles dihedrals

writepsf "my_1f7a.psf"
writepdb "my_1f7a.pdb"

foreach f $LOCALFILES {
  exec rm $f
}

mol delete top
mol new my_1f7a.psf
mol addfile my_1f7a.pdb
set a [atomselect top all]
set fix [atomselect top "noh"]
$a set beta 0
$fix set beta 1
$a writepdb my_1f7a_fix.pdb

exit
