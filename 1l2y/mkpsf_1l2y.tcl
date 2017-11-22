# VMD/psfgen script for generating psf/pdb pair for PDB 1l2y
# Trp-Cage engineered mini-protein
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

mol new 1l2y.pdb

set a [atomselect top "all"]
$a frame 0
$a writepdb "tmp-a.pdb"
lappend LOCALFILES tmp-a.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $PSFGEN_BASEDIR/charmm/prod.top
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str

pdbalias atom ILE CD1 CD
pdbalias residue HIS HSD

segment A {
  pdb tmp-a.pdb
}
coordpdb tmp-a.pdb A

guesscoord

regenerate angles dihedrals

writepsf "my_1l2y.psf"
writepdb "my_1l2y.pdb"

foreach f $LOCALFILES {
  exec rm $f
}

mol delete top
mol new my_1l2y.psf
mol addfile my_1l2y.pdb
set a [atomselect top all]
set fix [atomselect top "noh"]
$a set beta 0
$fix set beta 1
$a writepdb my_1l2y_fix.pdb

exit

