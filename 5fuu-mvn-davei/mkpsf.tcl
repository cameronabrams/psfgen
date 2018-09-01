# VMD/psfgen script for generated psf/pdb pair for a 5FUU-MVNDAVEI complex
#
# using a 5fuu from psfgen/5fuu and an MVNDAVEI from psfgen/2hhy
#
# cameron f abrams (c) 2013-2018
# drexel university
# chemical and biological engineering

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

source ${PSFGEN_BASEDIR}/src/loopmc.tcl

# a MVN* DAVEI with bound mannobiose (chain M) from 2hyy
mol new my_2yhh.psf
mol addfile my_2yhh_vac.coor
set mvn_id [molinfo top get id]
set a [atomselect $mvn_id "all"]
set p [atomselect $mvn_id "chain A"]
set m [atomselect $mvn_id "chain M and not name O1 HO1"]
$p set chain X
$p set segname X

# resids of the mannobiose you want to dock the MVN*DAVEI to
set mbs { 2452 2458 }
# 5fuu Env trimer, elaborated with MPER's and Man9's
# resids named in $mbs must exist in one of the Man9's
mol new my_5fuu.psf
mol addfile my_5fuu_vac.coor
set env_id [molinfo top get id]
set b [atomselect $env_id all]
$b writepdb t.pdb

package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf

set targ [atomselect $env_id "chain A and resid $mbs"]
$a move [measure fit $m $targ]
$targ set {x y z} [$m get {x y z}]
$p writepdb X.pdb

resetpsf
readpsf my_5fuu.psf
coordpdb t.pdb

segment X {
   pdb X.pdb
}
coordpdb X.pdb X
patch DISU X:8   X:24
patch DISU X:60  X:80
patch DISU X:63  X:78

writepsf "my_complex.psf"
writepdb "my_complex.pdb"

mol new my_complex.psf
mol addfile my_complex.pdb

set sel [atomselect top "chain X and resid 123 to 131"]
fold_alpha_helix top $sel extrabonds-X.inp

set a [atomselect top all]
$a set beta 0.0

set b [atomselect top "chain A B C D E F and not glycan and noh"]
$b set beta 1.0

$a writepdb "my_complex_fix.pdb"

# generate a colvars input file for vacuum MD SMD (vac stage1)
# the objective is to place the trp3 peptide part of DAVEI
# as close as possible to the nearest MPER segment

# Stage 1: drag Trp3's toward Env MPER's

set env_targ_as "664 to 672"
set env_targ_cv "664-672"
set env_avoid_cv "610-630"
set tp [lindex [$sel get {x y z}] 0]
set dx [veclength [vecsub $tp [measure center [atomselect top "chain F and resid $env_targ_as"]]]]
set fp [open "drag_daveis_colvars.inp" "w"]
puts $fp "
colvarstrajfrequency 100
scriptedcolvarforces on
colvar {
   name trp3_1
   distance {
     group1 {
       psfsegid X
       atomnameresiduerange CA 123-131
     }
     group2 {
       psfsegid F
       atomnameresiduerange CA $env_targ_cv
     }
   }
}
colvar {
   name trp3_1av
   distance {
     group1 {
       psfsegid X
       atomnameresiduerange CA 123-131
     }
     group2 {
       psfsegid D
       atomnameresiduerange CA $env_avoid_cv
     }
   }
}
harmonic {
  colvars trp3_1
  centers $dx
  targetcenters 4.0
  forceconstant 10.0
  targetnumsteps 60000
}
"
close $fp
quit


