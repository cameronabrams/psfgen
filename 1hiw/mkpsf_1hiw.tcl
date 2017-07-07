# VMD/psfgen script for generating psf/pdb pair for PDB 1hiw
# trimeric HIV-1 matrix protein (chains A, B, C)
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

mol new 1hiw.pdb
foreach c { A B C } {
  set a [atomselect top "protein and chain $c"]
  $a writepdb "${c}.pdb"; lappend LOCALFILES ${c}.pdb
  set w [atomselect top "water and same residue as within 4.0 of (protein and chain $c)"]
  $w set resname TIP3
  set wo [atomselect top "(water and same residue as within 4.0 of (protein and chain $c)) and name O"]
  $wo set name OH2
  $w writepdb "${c}W.pdb"; lappend LOCALFILES ${c}W.pdb
}

mol delete top

package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str

alias residue HIS HSD
alias atom ILE CD1 CD

foreach c { A B C } {
  segment $c {
    pdb ${c}.pdb
  }
  segment ${c}W {
    auto none
    pdb ${c}W.pdb
  }
  coordpdb ${c}.pdb $c
  coordpdb ${c}W.pdb ${c}W
}

guesscoord

writepsf "my_1hiw.psf"
writepdb "my_1hiw.pdb"

foreach f $LOCALFILES {
  exec rm $f
}

mol delete top
mol new my_1hiw.psf
mol addfile my_1hiw.pdb
set a [atomselect top all]
set fix [atomselect top "noh"]
$a set beta 0
$fix set beta 1
$a writepdb my_1hiw_fix.pdb

exit

