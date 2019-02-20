# VMD/psfgen script for generating psf/pdb pair for PDB 2k7w
# human BAX with BIM SAHB BH3-like peptide bound
#
# cameron f abrams (c) 2019
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

set BIM 0  ; # set to 1 to include the BIM SAHB peptide
set FRM 0  ; # allow user to select which of the frames in the NMR structure to use
set P168G 0 ; # allow user to perform P168G mutation
for { set a 0 } { $a < [llength $argv] } { incr a } {
   set arg [lindex $argv $a]
   if { $arg == "-bim" } {
      set BIM 1
   }
   if { $arg == "-frm" } {
      incr a
      set FRM [lindex $argv $a]
   }
   if { $arg == "-P168G" } {
      set P168G 1
   }
}

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

mol new 2k7w.pdb

set all [atomselect top "all"]
$all frame $FRM

set A [atomselect top "protein and chain A"]
$A writepdb A.pdb; lappend LOCALFILES A.pdb
if { $BIM == "1" } {
  set B [atomselect top "protein and chain B"]
  $B writepdb B.pdb; lappend LOCALFILES B.pdb
}

package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $PSFGEN_BASEDIR/charmm/prod.top

pdbalias atom ILE CD1 CD
pdbalias residue HIS HSD

segment A { 
  pdb A.pdb
  if { $P168G == "1" } {
     mutate 168 GLY
  }
}
coordpdb A.pdb A

if { $BIM == "1" } {
   segment B {
     pdb B.pdb
   }
   coordpdb B.pdb B
}

guesscoord

regenerate angles dihedrals

writepsf "my_2k7w.psf"
writepdb "my_2k7w.pdb"

foreach f $LOCALFILES {
  exec rm $f
}

mol delete top
mol new my_2k7w.psf
mol addfile my_2k7w.pdb
set a [atomselect top all]
set fix [atomselect top "noh"]
$a set beta 0
$fix set beta 1
$a writepdb my_2k7w_fix.pdb

exit
