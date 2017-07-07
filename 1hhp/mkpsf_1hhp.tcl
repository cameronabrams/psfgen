# VMD/psfgen script for generating psf/pdb pair for PDB 1hhp
# dimeric apo HIV-1 protease
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
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

mol new 1hhp.pdb

set a [atomselect top "all"]
$a set chain B
$a set segid B
$a set segname B
$a move { { 0.000000  1.000000  0.000000        0.00000} 
          { 1.000000  0.000000  0.000000        0.00000} 
          { 0.000000  0.000000  -1.000000       0.00000} 
          { 0 0 0 1 } }
$a writepdb B.pdb
lappend LOCALFILES B.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $PSFGEN_BASEDIR/charmm/prod.top
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str

pdbalias atom ILE CD1 CD
pdbalias residue HIS HSD

segment A {
  pdb 1hhp.pdb
}
coordpdb 1hhp.pdb A

segment B {
  pdb B.pdb
}
coordpdb B.pdb B

patch ASPP A:25
patch ASPP B:25
patch PROD A:44
patch PROD B:44

guesscoord

regenerate angles dihedrals

writepsf "my_1hhp.psf"
writepdb "my_1hhp.pdb"

foreach f $LOCALFILES {
  exec rm $f
}

mol delete top
mol new my_1hhp.psf
mol addfile my_1hhp.pdb
set a [atomselect top all]
set fix [atomselect top "noh"]
$a set beta 0
$fix set beta 1
$a writepdb my_1hhp_fix.pdb

exit

