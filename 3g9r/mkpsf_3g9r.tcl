# VMD/psfgen script for generating psf/pdb pair for PDB 3g9r
# trimeric construct of HIV-1 MPER
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

mol new 3g9r.pdb

set wat [atomselect top "water and within 3.0 of (protein and chain A B C and resid < 23)"]
$wat set name OH2
$wat set resname TIP3
$wat writepdb wat.pdb ; lappend LOCALFILES wat.pdb

set pro [atomselect top "protein and chain A and resid < 23"]
$pro writepdb pro_A.pdb; lappend LOCALFILES pro_A.pdb
set pro [atomselect top "protein and chain B and resid < 23"]
$pro writepdb pro_B.pdb; lappend LOCALFILES pro_B.pdb
set pro [atomselect top "protein and chain C and resid < 23"]
$pro writepdb pro_C.pdb; lappend LOCALFILES pro_C.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str

pdbalias atom ILE CD1 CD
pdbalias residue HIS HSD

segment A {
  pdb pro_A.pdb
  first NNEU
  last CNEU
}
segment B {
  pdb pro_B.pdb
  first NNEU
  last CNEU
}
segment C {
  pdb pro_C.pdb
  first NNEU
  last CNEU
}
segment WX {
  auto none
  pdb wat.pdb
}

coordpdb pro_A.pdb A
coordpdb pro_B.pdb B
coordpdb pro_C.pdb C
coordpdb wat.pdb WX

guesscoord

regenerate angles dihedrals

writepsf "my_3g9r.psf"
writepdb "my_3g9r.pdb"

foreach f $LOCALFILES {
  exec rm $f
}

mol delete top
mol new my_3g9r.psf
mol addfile my_3g9r.pdb
set a [atomselect top all]
set fix [atomselect top "noh"]
$a set beta 0
$fix set beta 1
$a writepdb my_3g9r_fix.pdb

exit

