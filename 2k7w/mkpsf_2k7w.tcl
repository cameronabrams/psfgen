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

set BIM 1  ; # set to 0 to not include the BIM SAHB peptide
set FRM 0  ; # allow user to select which of the frames in the NMR structure to use
set P168G 0 ; # allow user to perform P168G mutation
set cisP168 0 ; # allow user to make the 167-168 peptide bond cis
for { set a 0 } { $a < [llength $argv] } { incr a } {
   set arg [lindex $argv $a]
   if { $arg == "-nobim" } {
      set BIM 0
   }
   if { $arg == "-frm" } {
      incr a
      set FRM [lindex $argv $a]
   }
   if { $arg == "-P168G" } {
      set P168G 1
   }
   if { $arg == "-cisP168" } {
      set cisP168 1
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

if { $cisP168 == "1" } {
   set fp [open "go-stage-2" "w"]
   puts $fp "go for cis!"
   close $fp
   set CA1 [atomselect top "protein and chain A and resid 167 and name CA"]
   set C   [atomselect top "protein and chain A and resid 167 and name C"]
   set N   [atomselect top "protein and chain A and resid 168 and name N"]
   set CA2 [atomselect top "protein and chain A and resid 168 and name CA"]
   set fp [open "cvomega.inp" "w"]
   puts $fp "
# CV input for trans-cis isomerization of 167-168 peptide bond
colvarsTrajFrequency 100

colvar {
  name omega
  dihedral {
     group1 {
        atomNumbers [$CA1 get serial]
     }
     group2 {
        atomNumbers [$C get serial]
     }
     group3 {
        atomNumbers [$N get serial]
     }
     group4 {
        atomNumbers [$CA2 get serial]
     }
  } 
}

harmonic {
  name make_cis
  colvars omega
  forceConstant 10.0
  centers 180.0
  targetCenters 0.0
  targetNumSteps 20000
}
"
close $fp
}

exit
