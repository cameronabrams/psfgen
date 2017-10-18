# VMD/psfgen script for generating psf/pdb pair for PDB 5ggr
# nivolumab ("opdivo"/bms) in complex with PD-1
# 
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

set DOMC 1

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

mol new 5ggr.pdb

# extract contiguous protein segments as individual pdb's                                                                                                             
set segs { 
  { H 2 128 } { H 134 213 } { L 1 212 } { Z 27 71 } { Z 75 146 }
}

set gaps { 
  { H 129 133 } { Z 72 74 }
}

set ns {}
set mln { 
  { H 129 128 } { Z 72 71 }
}
foreach m $mln {
  lappend ns [cacoIn_nOut [lindex $m 2] [lindex $m 0] 0]
}
foreach s $segs {
  [atomselect top "protein and chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] writepdb "[lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb"
  lappend LOCALFILES [lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb
}

mol delete top
package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

segment H {
  pdb H_2_to_128.pdb
  residue 129 LYS H
  residue 130 SER H
  residue 131 THR H
  residue 132 SER H
  residue 133 GLY H  
  pdb H_134_to_213.pdb
}

segment L {
  pdb L_1_to_212.pdb
}

segment Z {
  pdb Z_27_to_71.pdb
  residue 72 PRO Z
  residue 73 SER Z
  residue 74 ASN Z
  pdb Z_75_to_146.pdb
}

foreach s $segs {
  coordpdb [lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb [lindex $s 0]
}

# set N positions
foreach m $mln n $ns {
  puts "coord [lindex $m 0] [lindex $m 1] N $n"
  coord [lindex $m 0] [lindex $m 1] N $n
}

# disulfides
patch DISU H:22 H:96
patch DISU H:140 H:196
patch DISU L:23 L:88
patch DISU Z:54 Z:123

guesscoord
regenerate angles dihedrals

writepsf "my_5ggr.psf"
writepdb "unrelaxed.pdb"

lappend LOCALFILES unrelaxed.pdb
package require Orient
namespace import Orient::orient

mol delete top
mol new my_5ggr.psf
mol addfile unrelaxed.pdb
set molid [molinfo top get id]
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
set fab [atomselect top "chain H L"]
set I [Orient::calc_principalaxes $fab]
set A [Orient::orient $fab [lindex $I 2] {0 0 -1}]
$a move $A
$a writepdb "unrelaxed2.pdb"
lappend LOCALFILES "unrelaxed2.pdb"
 
if { $DOMC == "1" } {
 set nc 1000
 set rcut 3.0
 set temperature 2.5
 set k 10.0
 set r0 1.5
 set bg [atomselect ${molid} "noh"]
 foreach l $gaps {
   set chain [lindex $l 0]
   set residueList [[atomselect ${molid} "chain $chain and resid [lindex $l 1] to [lindex $l 2] and name CA"] get residue]
   do_loop_mc ${residueList} ${chain} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999]
 }
}                                                                                                                                                                     
$a writepdb "my_5ggr_mcOut.pdb"

# make a pdb file that fixes all heavy atoms in the original
# crystal structure -- all added atoms are set as unfixed
# for a minimization
mol delete top
mol new my_5ggr.psf
mol addfile unrelaxed2.pdb
set a [atomselect top all]
$a set beta 0
foreach s $segs {
  [atomselect top "chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] set beta 1
}
[atomselect top "not noh"] set beta 0

$a writepdb "my_5ggr_fix.pdb"

# clean up
foreach f $LOCALFILES {
  exec rm $f
}

quit
