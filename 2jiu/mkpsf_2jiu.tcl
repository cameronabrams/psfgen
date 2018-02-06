# VMD/psfgen script for generating psf/pdb pair for PDB 2jiu
# EGFRK with ATP-like ligand bound;
# Note, we replace this ATP-like ligand with MgATP
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

# turn on MC loop conformation sampling
set DOMC 1

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

mol new 2jiu.pdb
mol new 2gs6.pdb

# define binding site based on the ADP-bound 2gs6 structure; numbering
# is different in 2jiu
set gs6_bs [atomselect 1 "chain A and backbone and resid 721 768 831"]
set jiu_bs [atomselect 0 "chain A and backbone and resid 745 792 855"]
# select that part of the ATP-like analog that is common to ATP
# and fit it by binding-site alignment into the 2jiu coordinate system
# and save it
set analog [atomselect 1 "resname 112 and not (name S1G C1S C2S O2S NS)"]
$analog move [measure fit $gs6_bs $jiu_bs]
$analog set resname "ATP"
$analog set chain "Q"
$analog set resid 1
$analog writepdb 2gs6_resname112_aligned.pdb
mol delete 0
mol delete 1
lappend LOCALFILES 2gs6_resname112_aligned.pdb

# reload; perform a two-stage alignment of the A-chain ATP molecule
# from 5uds to the ATP fragment already in the 2jiu coordinate frame
# first stage aligns the adenosine; second aligns the Mg2+ into
# the pseudo six-membered ring formed by O2B O3B O3G PB PG.
# (a single stage did not work robustly in my testing)
mol new 2gs6_resname112_aligned.pdb
set id1 [molinfo top get id]
set a [atomselect $id1 "name C4' O4' C1' C5 N7 C8 N9 N1 C2 N3 C4 C6 N6 C2' O2'"] 
mol new 5uds.pdb
set id2 [molinfo top get id]
set atpmg [atomselect $id2 "chain A and resname ATP MG"]
set b [atomselect $id2 "chain A and resname ATP and name C4' O4' C1' C5 N7 C8 N9 N1 C2 N3 C4 C6 N6 C2' O2'"]
$atpmg move [measure fit $b $a]
set a [atomselect $id1 "name O2B O3B O3G PB PG"]
set b [atomselect $id2 "chain A and name O2B O3B O3G PB PG"]
set m [atomselect $id2 "chain A and resname MG"]
$m move [measure fit $b $a]
$m writepdb "MG.pdb"
lappend LOCALFILES MG.pdb
mol delete $id1
mol delete $id2

lappend LOCALFILE MG.pdb

mol new 2jiu.pdb
set molid [molinfo top get id]

set segs { { A 696 992 } { A 1010 1016 } }
set loops { { A 993 1009 } }
set ns {}
set mln { { A 993 992 } }

foreach m $mln {
   lappend ns [cacoIn_nOut [lindex $m 2] [lindex $m 0] $molid]
}

foreach s $segs {
  [atomselect $molid "protein and chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] writepdb "[lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb"
  lappend LOCALFILES [lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb
}

set wat [atomselect $molid "chain A and water"]
$wat set name OH2
$wat set resname TIP3
$wat set segname WX
$wat set chain WX
$wat writepdb "A_water.pdb"
lappend LOCALFILES A_water.pdb

mol delete $molid

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str
topology $env(HOME)/charmm/toppar/stream/na/toppar_all36_na_nad_ppi.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

segment Q {
  pdb 2gs6_resname112_aligned.pdb
}

segment MG {
  pdb MG.pdb
}

coordpdb 2gs6_resname112_aligned.pdb Q
coordpdb MG.pdb MG

segment A {
  pdb A_696_to_992.pdb
  residue 993 THR A
  residue 994 ASP A
  residue 995 SER A
  residue 996 ASN A
  residue 997 PHE A
  residue 998 TYR A
  residue 999 ARG A
  residue 1000 ALA A
  residue 1001 LEU A
  residue 1002 MET A
  residue 1003 ASP A
  residue 1004 GLU A
  residue 1005 GLU A
  residue 1006 ASP A
  residue 1007 MET A
  residue 1008 ASP A
  residue 1009 ASP A
  pdb A_1010_to_1016.pdb
}

segment WX {
  auto none
  pdb A_water.pdb
}

coordpdb A_696_to_992.pdb A
coordpdb A_1010_to_1016.pdb A
coordpdb A_water.pdb WX

foreach m $mln n $ns {
  puts "coord [lindex $m 0] [lindex $m 1] N $n"
  coord [lindex $m 0] [lindex $m 1] N $n
}

guesscoord

regenerate angles dihedrals

writepsf "my_2jiu.psf"
writepdb "unrelaxed.pdb"

lappend LOCALFILES unrelaxed.pdb

mol new my_2jiu.psf
mol addfile unrelaxed.pdb
set molid [molinfo top get id]
set a [atomselect $molid all]
set or [measure center $a weight mass]
$a moveby [vecscale -1 $or]

if { $DOMC == "1" } {
  set nc 1000
  set rcut 5.0
  set temperature 1.0
  set k 1.0
  set r0 1.5
  set bg [atomselect ${molid} "noh"]
  foreach l $loops {
    set chain [lindex $l 0]
    set residueList [[atomselect ${molid} "chain $chain and resid [lindex $l 1] to [lindex $l 2] and name CA"] get residue]
    do_loop_mc ${residueList} ${chain} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] 
  }
}

$a writepdb "my_2jiu_mcOut.pdb"
mol delete $molid
mol new my_2jiu.psf
mol addfile my_2jiu_mcOut.pdb
set a [atomselect top all]
$a set beta 0
foreach s $segs {
  [atomselect top "chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] set beta 1
}
[atomselect top "name OH2"] set beta 1
[atomselect top "chain A and resid 992"] set beta 0
[atomselect top "not noh"] set beta 0
[atomselect top "resname ATP MG"] set beta 0

$a writepdb "my_2jiu_fix.pdb"

# clean up
foreach f $LOCALFILES { 
  exec rm $f
}

exit

