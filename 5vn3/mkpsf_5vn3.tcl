# VMD/psfgen script for generating psf/pdb pair for PDB 5vn3
# trimeric HIV-1 gp140 in "open" state
# 
# note: 17b antibodies (6 chains) and sCD4's (3 chains) are 
# not included as yet.
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

mol new 5vn3.pdb

# extract contiguous protein segments as individual pdb's                                                                                                             
set segs { 
  { G 32 131 } { G 193 309 } { G 325 404 } { G 412 505 } { A 515 548 } { A 563 664}  
  { I 32 131 } { I 193 309 } { I 325 404 } { I 412 505 } { B 515 548 } { B 563 664} 
  { J 32 131 } { J 193 309 } { J 325 404 } { J 412 505 } { D 515 548 } { D 563 664} 
}

set gaps { 
  { G 132 192 } { G 311 324 } { G 407 411 } { A 549 562 }
  { I 132 192 } { I 311 324 } { I 407 411 } { B 549 562 }
  { J 132 192 } { J 311 324 } { J 407 411 } { D 549 562 }
}

set ns {}
set mln { 
  { G 132 131 } { G 311 309 } { G 407 404 } { A 549 548 }
  { I 132 131 } { I 311 309 } { I 407 404 } { B 549 548 }
  { J 132 131 } { J 311 309 } { J 407 404 } { D 549 548 }
}
foreach m $mln {
  lappend ns [cacoIn_nOut [lindex $m 2] [lindex $m 0] 0]
}
foreach s $segs {
  [atomselect top "protein and chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] writepdb "[lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb"
  lappend LOCALFILES [lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb
}

# save each chain's glycan residues with CHARMM-compatible name
foreach chain { G I J A B D } {
  set g [atomselect top "chain $chain and not protein"]
  [atomselect top "chain $chain and resname NAG"] set resname BGNA
  [atomselect top "chain $chain and resname MAN"] set resname AMAN
  [atomselect top "chain $chain and resname BMA"] set resname BMAN
  [atomselect top "chain $chain and resname FUC"] set resname AFUC
  [atomselect top "chain $chain and resname GAL"] set resname BGAL
  $g writepdb "${chain}S.pdb"
  lappend LOCALFILES ${chain}S.pdb
}


mol delete top
package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias atom BGNA N2 N
pdbalias atom BGNA C7 C
pdbalias atom BGNA O7 O
pdbalias atom BGNA C8 CT

segment G {
  pdb G_32_to_131.pdb
  residue 132 GLY G
  residue 133 GLY G
  residue 134 GLY G; # this is just because I don't want to try to model in 63 missing residues of V1/2!
  pdb G_193_to_309.pdb
  residue 311 GLY G
  residue 312 PRO G
  residue 313 GLY G
  residue 314 ARG G
  residue 315 ALA G
  residue 316 PHE G
  residue 317 TYR G
  residue 318 ALA G
  residue 319 THR G
  residue 320 GLY G
  residue 321 ASP G
  residue 322 ILE G
  residue 323 ILE G
  residue 324 GLY G
  pdb G_325_to_404.pdb
  residue 407 PRO G
  residue 408 THR G
  residue 409 GLY G
  residue 410 GLY G
  residue 411 GLU G
  pdb G_412_to_505.pdb
}

segment GS {
  pdb GS.pdb
}

segment A {
  pdb A_515_to_548.pdb
  residue 549 VAL A
  residue 550 GLN A
  residue 551 GLN A
  residue 552 GLN A
  residue 553 ASN A
  residue 554 ASN A
  residue 555 LEU A
  residue 556 LEU A
  residue 557 ARG A
  residue 558 ALA A
  residue 559 PRO A
  residue 560 GLU A
  residue 561 ALA A
  residue 562 GLN A
  pdb A_563_to_664.pdb
}

segment AS {
  pdb AS.pdb
}

segment I {
  pdb I_32_to_131.pdb
  residue 132 GLY I
  residue 133 GLY I
  residue 134 GLY I; # this is just because I don't want to try to model in 63 missing residues of V1/2!
  pdb I_193_to_309.pdb
  residue 311 GLY I
  residue 312 PRO I
  residue 313 GLY I
  residue 314 ARG I
  residue 315 ALA I
  residue 316 PHE I
  residue 317 TYR I
  residue 318 ALA I
  residue 319 THR I
  residue 320 GLY I
  residue 321 ASP I
  residue 322 ILE I
  residue 323 ILE I
  residue 324 GLY I
  pdb I_325_to_404.pdb
  residue 407 PRO I
  residue 408 THR I
  residue 409 GLY I
  residue 410 GLY I
  residue 411 GLU I
  pdb I_412_to_505.pdb
}

segment IS {
  pdb IS.pdb
}

segment B {
  pdb B_515_to_548.pdb
  residue 549 VAL B
  residue 550 GLN B
  residue 551 GLN B
  residue 552 GLN B
  residue 553 ASN B
  residue 554 ASN B
  residue 555 LEU B
  residue 556 LEU B
  residue 557 ARG B
  residue 558 ALA B
  residue 559 PRO B
  residue 560 GLU B
  residue 561 ALA B
  residue 562 GLN B
  pdb B_563_to_664.pdb
}

segment BS {
  pdb BS.pdb
}

segment J {
  pdb J_32_to_131.pdb
  residue 132 GLY J
  residue 133 GLY J
  residue 134 GLY J; # this is just because I don't want to try to model in 63 missing residues of V1/2!
  pdb J_193_to_309.pdb
  residue 311 GLY J
  residue 312 PRO J
  residue 313 GLY J
  residue 314 ARG J
  residue 315 ALA J
  residue 316 PHE J
  residue 317 TYR J
  residue 318 ALA J
  residue 319 THR J
  residue 320 GLY J
  residue 321 ASP J
  residue 322 ILE J
  residue 323 ILE J
  residue 324 GLY J
  pdb J_325_to_404.pdb
  residue 407 PRO J
  residue 408 THR J
  residue 409 GLY J
  residue 410 GLY J
  residue 411 GLU J
  pdb J_412_to_505.pdb
}

segment JS {
  pdb JS.pdb
}

segment D {
  pdb D_515_to_548.pdb
  residue 549 VAL D
  residue 550 GLN D
  residue 551 GLN D
  residue 552 GLN D
  residue 553 ASN D
  residue 554 ASN D
  residue 555 LEU D
  residue 556 LEU D
  residue 557 ARG D
  residue 558 ALA D
  residue 559 PRO D
  residue 560 GLU D
  residue 561 ALA D
  residue 562 GLN D
  pdb D_563_to_664.pdb
}

segment DS {
  pdb DS.pdb
}

foreach s $segs {
  coordpdb [lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb [lindex $s 0]
}

# coordpdb sugars
coordpdb GS.pdb GS
coordpdb IS.pdb IS
coordpdb JS.pdb JS
coordpdb AS.pdb AS
coordpdb BS.pdb BS
coordpdb DS.pdb DS

# set N positions
foreach m $mln n $ns {
  puts "coord [lindex $m 0] [lindex $m 1] N $n"
  coord [lindex $m 0] [lindex $m 1] N $n
}

# disulfides
patch DISU G:54 G:74
patch DISU G:119 G:205
patch DISU G:126 G:196
# patch DISU G:131 G:157
patch DISU G:218 G:247
patch DISU G:228 G:239
patch DISU G:296 G:331
patch DISU G:378 G:445
patch DISU G:385 G:418
patch DISU G:501 A:605; # sosip!
patch DISU A:598 A:604
patch DISU I:54 I:74
patch DISU I:119 I:205
patch DISU I:126 I:196
# patch DISU I:131 I:157
patch DISU I:218 I:247
patch DISU I:228 I:239
patch DISU I:296 I:331
patch DISU I:378 I:445
patch DISU I:385 I:418
patch DISU I:501 B:605; # sosip!
patch DISU B:598 B:604
patch DISU J:54 J:74
patch DISU J:119 J:205
patch DISU J:126 J:196
# patch DISU J:131 J:157
patch DISU J:218 J:247
patch DISU J:228 J:239
patch DISU J:296 J:331
patch DISU J:378 J:445
patch DISU J:385 J:418
patch DISU J:501 D:605; # sosip!
patch DISU D:598 D:604

patch NGLB A:611 AS:701
patch NGLB A:616 AS:702
patch NGLB A:625 AS:703
patch NGLB A:637 AS:704
patch NGLB B:611 BS:701
patch NGLB B:616 BS:702
patch NGLB B:625 BS:703
patch NGLB B:637 BS:704
patch NGLB D:611 DS:701
patch NGLB D:616 DS:702
patch NGLB D:625 DS:703
patch NGLB D:637 DS:704

patch NGLB G:88 GS:601
patch NGLB G:234 GS:602
patch 14bb GS:602 GS:603
patch NGLB G:241 GS:604
patch 14bb GS:604 GS:605
patch NGLB G:262 GS:607
patch 14bb GS:607 GS:608
patch 14bb GS:608 GS:609
patch 13ab GS:609 GS:610
patch 16ab GS:609 GS:613
patch 12ab GS:610 GS:611
patch 12ab GS:611 GS:612
patch NGLB G:276 GS:614
patch 14bb GS:614 GS:615
patch 14bb GS:615 GS:616
patch 16bb GS:616 GS:617
patch NGLB G:295 GS:618
patch 14bb GS:618 GS:619
patch NGLB G:332 GS:620
patch 14bb GS:620 GS:621
patch NGLB G:339 GS:622
patch NGLB G:355 GS:623
patch NGLB G:362 GS:624
patch 14bb GS:624 GS:625
patch 14bb GS:625 GS:626
patch 13bb GS:626 GS:627
patch NGLB G:386 GS:628
patch 14bb GS:628 GS:629
patch 14bb GS:629 GS:630
patch NGLB G:392 GS:637
patch 14bb GS:637 GS:638
patch 14bb GS:638 GS:639
patch NGLB G:397 GS:631
patch NGLB G:413 GS:632
patch 14bb GS:632 GS:633
patch 14bb GS:633 GS:634
patch NGLB G:448 GS:635
patch 14bb GS:635 GS:636

patch NGLB I:88 IS:601
patch NGLB I:234 IS:602
patch 14bb IS:602 IS:603
patch NGLB I:241 IS:604
patch 14bb IS:604 IS:605
patch NGLB I:262 IS:607
patch 14bb IS:607 IS:608
patch 14bb IS:608 IS:609
patch 13ab IS:609 IS:610
patch 16ab IS:609 IS:613
patch 12ab IS:610 IS:611
patch 12ab IS:611 IS:612
patch NGLB I:276 IS:614
patch 14bb IS:614 IS:615
patch 14bb IS:615 IS:616
patch 16bb IS:616 IS:617
patch NGLB I:295 IS:618
patch 14bb IS:618 IS:619
patch NGLB I:332 IS:620
patch 14bb IS:620 IS:621
patch NGLB I:339 IS:622
patch NGLB I:355 IS:623
patch NGLB I:362 IS:624
patch 14bb IS:624 IS:625
patch 14bb IS:625 IS:626
patch 13bb IS:626 IS:627
patch NGLB I:386 IS:628
patch 14bb IS:628 IS:629
patch 14bb IS:629 IS:630
patch NGLB I:392 IS:637
patch 14bb IS:637 IS:638
patch 14bb IS:638 IS:639
patch NGLB I:397 IS:631
patch NGLB I:413 IS:632
patch 14bb IS:632 IS:633
patch 14bb IS:633 IS:634
patch NGLB I:448 IS:635
patch 14bb IS:635 IS:636

patch NGLB J:88 JS:601
patch NGLB J:234 JS:602
patch 14bb JS:602 JS:603
patch NGLB J:241 JS:604
patch 14bb JS:604 JS:605
patch NGLB J:262 JS:607
patch 14bb JS:607 JS:608
patch 14bb JS:608 JS:609
patch 13ab JS:609 JS:610
patch 16ab JS:609 JS:613
patch 12ab JS:610 JS:611
patch 12ab JS:611 JS:612
patch NGLB J:276 JS:614
patch 14bb JS:614 JS:615
patch 14bb JS:615 JS:616
patch 16bb JS:616 JS:617
patch NGLB J:295 JS:618
patch 14bb JS:618 JS:619
patch NGLB J:332 JS:620
patch 14bb JS:620 JS:621
patch NGLB J:339 JS:622
patch NGLB J:355 JS:623
patch NGLB J:362 JS:624
patch 14bb JS:624 JS:625
patch 14bb JS:625 JS:626
patch 13bb JS:626 JS:627
patch NGLB J:386 JS:628
patch 14bb JS:628 JS:629
patch 14bb JS:629 JS:630
patch NGLB J:392 JS:637
patch 14bb JS:637 JS:638
patch 14bb JS:638 JS:639
patch NGLB J:397 JS:631
patch NGLB J:413 JS:632
patch 14bb JS:632 JS:633
patch 14bb JS:633 JS:634
patch NGLB J:448 JS:635
patch 14bb JS:635 JS:636

guesscoord
regenerate angles dihedrals

writepsf "my_5vn3.psf"
writepdb "unrelaxed.pdb"

lappend LOCALFILES unrelaxed.pdb

mol delete top
mol new my_5vn3.psf
mol addfile unrelaxed.pdb
set molid [molinfo top get id]
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
set ca [measure center [atomselect top "protein and chain G I J"] weight mass]
set cb [measure center [atomselect top "protein and chain A B D"] weight mass]
set pi 3.1415928
set dv [vecsub $ca $cb]
set d [veclength $dv]
set cp [expr [lindex $dv 0]/$d]
set sp [expr [lindex $dv 1]/$d]
set p [expr acos($cp)]
if {[expr $sp < 0.0]} {
  set p [expr 2*$pi-$p]
}
set ct [expr [lindex $dv 2]/$d]
set t [expr acos($ct)]
$a move [transaxis z [expr -1 * $p] rad]
$a move [transaxis y [expr -1 * $t] rad]
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
$a writepdb "my_5vn3_mcOut.pdb"

# make a pdb file that fixes all heavy atoms in the original
# crystal structure -- all added atoms are set as unfixed
# for a minimization
mol delete top
mol new my_5vn3.ps
mol addfile unrelaxed2.pdb
set a [atomselect top all]
$a set beta 0
foreach s $segs {
  [atomselect top "chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] set beta 1
}
[atomselect top "not noh"] set beta 0

$a writepdb "my_5vn3_fix.pdb"

# clean up
foreach f $LOCALFILES {
  exec rm $f
}

quit


resetpsf

# make a fixed-atoms pdb for relaxation
# mol new my_4zmj_protomer_GB_glycans.psf
# mol addfile my_4zmj_protomer_GB_glycans_rawloops_tmp1.pdb

# set a [atomselect top all]
#$a set beta 1 ;# pretty much everything is fixed except...

#[atomselect top "chain G and resid 185 to 187"] set beta 0
#[atomselect top "chain G and resid 397 to 409"] set beta 0
#[atomselect top "chain B and resid 512 to 521"] set beta 0
#[atomselect top "chain B and resid 547 to 569"] set beta 0
#$a writepdb "my_fix_GB.pdb"

#mol delete top

foreach g {G E F} b {B C D} {
  readpsf my_4zmj_protomer_${g}${b}_glycans.psf pdb my_4zmj_protomer_${g}${b}_glycans_rawloops_tmp1.pdb
}
writepsf "my_4zmj.psf"
writepdb "trimer.pdb"
lappend LOCALFILES trimer.pdb

mol new my_4zmj.psf
mol addfile trimer.pdb
set molid [molinfo top get id]

set a [atomselect top all]
$a set beta 1 ;# pretty much everything is fixed except...

[atomselect top "chain G E F and resid 185 to 187"] set beta 0
[atomselect top "chain G E F and resid 397 to 409"] set beta 0
[atomselect top "chain B C D and resid 512 to 521"] set beta 0
[atomselect top "chain B C D and resid 547 to 569"] set beta 0
$a writepdb "my_4zmj_fix.pdb"

$a set beta 0

foreach b { B C D } {
  set residueList [[atomselect top "chain ${b} and resid 547 to 568 and name CA"] get residue]
  set r1 [lindex $residueList 0]
  set r2 [lindex $residueList end]

  Crot_phi $r1 $r2 ${b} $molid -180
  Crot_psi $r1 $r2 ${b} $molid -30
  Crot_phi $r1 $r2 ${b} $molid 10
}

$a writepdb "trimer_mcIn.pdb"
lappend LOCALFILES trimer_mcIn.pdb

set nc 1000
set rcut 3.0
set temperature 2.5
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]
foreach g {G E F} b {B C D} {
  set residueList [[atomselect ${molid} "chain ${g} and resid 186 and name CA"] get residue]
  do_loop_mc ${residueList} ${g} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999]
  set residueList [[atomselect ${molid} "chain ${g} and resid 398 to 408 and name CA"] get residue]
  do_loop_mc ${residueList} ${g} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999]
  set residueList [[atomselect ${molid} "chain ${b} and resid 548 to 568 and name CA"] get residue]
  do_loop_mc ${residueList} ${b} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999]
}

$a writepdb "my_4zmj_mcOut.pdb"

foreach f $LOCALFILES {
  exec rm $f
}

quit

