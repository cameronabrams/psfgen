# VMD/psfgen script for generating psf/pdb pair for PDB 4zmj
# unliganded trimeric HIV-1 gp140
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
source ${PSFGEN_BASEDIR}/tcl/loopmc.tcl
set LOCALFILES {}

mol new 4zmj.pdb
set a [atomselect top all]
set ap [atomselect top "protein"]
set ag [atomselect top "not protein"]

$ap writepdb "GBp.pdb"; lappend LOCALFILES GBp.pdb
$ag writepdb "GBg.pdb"; lappend LOCALFILES GBg.pdb

set g [atomselect top "chain G"]
set b [atomselect top "chain B"]

set n186pos(G) [cacoIn_nOut 185 G 0]
set n398pos(G) [cacoIn_nOut 397 G 0]
set n548pos(B) [cacoIn_nOut 547 B 0]


set tmat3 {{-0.500000 -0.866025  0.000000      107.18000} {0.866025 -0.500000  0.000000      185.64121} {0.000000  0.000000  1.000000        0.00000} {0 0 0 1}}

$a move $tmat3
$g set chain E
$b set chain C
set n186pos(E) [cacoIn_nOut 185 E 0]
set n398pos(E) [cacoIn_nOut 397 E 0]
set n548pos(C) [cacoIn_nOut 547 C 0]
$ap writepdb "ECp.pdb"; lappend LOCALFILES ECp.pdb
$ag writepdb "ECg.pdb"; lappend LOCALFILES ECg.pdb

$a move $tmat3
$g set chain F
$b set chain D
set n186pos(F) [cacoIn_nOut 185 F 0]
set n398pos(F) [cacoIn_nOut 397 F 0]
set n548pos(D) [cacoIn_nOut 547 D 0]
$ap writepdb "FDp.pdb"; lappend LOCALFILES FDp.pdb
$ag writepdb "FDg.pdb"; lappend LOCALFILES FDg.pdb

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

foreach g {G E F} b {B C D} {
 
  mol new "${g}${b}p.pdb"
  set s1 [atomselect top "resid 34 to 185"]
  $s1 writepdb "${g}-s1.pdb"; lappend LOCALFILES ${g}-s1.pdb
  set s2 [atomselect top "resid 187 to 397"]
  $s2 writepdb "${g}-s2.pdb"; lappend LOCALFILES ${g}-s2.pdb
  set s3 [atomselect top "resid 411 to 505"]
  $s3 writepdb "${g}-s3.pdb"; lappend LOCALFILES ${g}-s3.pdb

  set s4 [atomselect top "resid 521 to 547"]
  $s4 writepdb "${b}-s4.pdb"; lappend LOCALFILES ${b}-s4.pdb
  set s5 [atomselect top "resid 569 to 664"]
  $s5 writepdb "${b}-s5.pdb"; lappend LOCALFILES ${b}-s5.pdb

  mol delete top
  mol new "${g}${b}g.pdb"
  [atomselect top "resname NAG"] set resname BGNA
  [atomselect top "resname MAN"] set resname AMAN
  [atomselect top "resname BMA"] set resname BMAN

  set ggg [atomselect top "chain $g"]
  $ggg writepdb "${g}g.pdb"; lappend LOCALFILES ${g}g.pdb
  set bbb [atomselect top "chain $b"]
  $bbb writepdb "${b}g.pdb"; lappend LOCALFILES ${b}g.pdb
  mol delete top

  resetpsf

  segment ${g} {
    pdb ${g}-s1.pdb
    residue 186A GLU ${g}
    residue 186B ASN ${g}
    residue 186C GLN ${g}
    residue 186D GLY ${g}
    residue 186E ASN ${g}
    residue 186F ARG ${g}
    residue 186G SER ${g}
    residue 186H ASN ${g}
    residue 186I ASN ${g}
    pdb ${g}-s2.pdb
    residue 398 THR ${g}
    residue 399 SER ${g}
    residue 400 VAL ${g}
    residue 401 GLN ${g}
    residue 402 GLY ${g}
    residue 403 SER ${g}
    residue 404 ASN ${g}
    residue 405 SER ${g}
    residue 406 THR ${g}
    residue 407 GLY ${g}
    residue 408 SER ${g}
    pdb ${g}-s3.pdb
  }

  segment ${g}S {
    pdb ${g}g.pdb
  }

  segment ${b} {
    residue 512 ALA ${b}
    residue 513 VAL ${b}
    residue 514 GLY ${b}
    residue 515 ILE ${b}
    residue 516 GLY ${b}
    residue 517 ALA ${b}
    residue 518 VAL ${b}
    residue 519 PHE ${b}
    residue 520 LEU ${b}
    pdb ${b}-s4.pdb
    residue 548 ILE ${b}
    residue 549 VAL ${b}
    residue 550 GLN ${b}
    residue 551 GLN ${b}
    residue 552 GLN ${b}
    residue 553 SER ${b}
    residue 554 ASN ${b}
    residue 555 LEU ${b}
    residue 556 LEU ${b}
    residue 557 ARG ${b}
    residue 558 ALA ${b}
    residue 559 PRO ${b}
    residue 560 GLU ${b}
    residue 561 ALA ${b}
    residue 562 GLN ${b}
    residue 563 GLN ${b}
    residue 564 HSD ${b}
    residue 565 LEU ${b}
    residue 566 LEU ${b}
    residue 567 LYS ${b}
    residue 568 LEU ${b}
    pdb ${b}-s5.pdb 
  }
  segment ${b}S {
    pdb ${b}g.pdb
  }

  coordpdb ${g}-s1.pdb ${g}
  coordpdb ${g}-s2.pdb ${g}
  coordpdb ${g}-s3.pdb ${g}
  coordpdb ${g}g.pdb ${g}S
  coordpdb ${b}-s4.pdb ${b}
  coordpdb ${b}-s5.pdb ${b}
  coordpdb ${b}g.pdb ${b}S

  coord ${g} 186A N $n186pos(${g})
  coord ${g} 398  N $n398pos(${g})
  coord ${b} 548  N $n548pos(${b})

  patch DISU ${g}:54 ${g}:74
  patch DISU ${g}:119 ${g}:205
  patch DISU ${g}:126 ${g}:196
  patch DISU ${g}:131 ${g}:157
  patch DISU ${g}:218 ${g}:247
  patch DISU ${g}:228 ${g}:239
  patch DISU ${g}:296 ${g}:331
  patch DISU ${g}:378 ${g}:445
  patch DISU ${g}:385 ${g}:418
  patch DISU ${g}:501 ${b}:605
  patch DISU ${b}:598 ${b}:604

  patch NGLB ${g}:156 ${g}S:615
  patch NGLB ${g}:160 ${g}S:616
  patch NGLB ${g}:197 ${g}S:617
  patch NGLB ${g}:234 ${g}S:601
  patch NGLB ${g}:262 ${g}S:602
  patch NGLB ${g}:276 ${g}S:608
  patch NGLB ${g}:295 ${g}S:619
  patch NGLB ${g}:301 ${g}S:620
  patch NGLB ${g}:332 ${g}S:621
  patch NGLB ${g}:339 ${g}S:609
  patch NGLB ${g}:355 ${g}S:610
  patch NGLB ${g}:363 ${g}S:611
  patch NGLB ${g}:386 ${g}S:612
  patch NGLB ${g}:392 ${g}S:613
  patch NGLB ${g}:448 ${g}S:614
  patch NGLB ${b}:611 ${b}S:701
  patch NGLB ${b}:618 ${b}S:702
  patch NGLB ${b}:637 ${b}S:703

  patch 14bb ${g}S:602 ${g}S:603
  patch 14bb ${g}S:617 ${g}S:618
  patch 14bb ${g}S:621 ${g}S:622

  patch 14bb ${g}S:603 ${g}S:604

  patch 13ab ${g}S:604 ${g}S:606
  patch 16ab ${g}S:604 ${g}S:605
  patch 12ab ${g}S:606 ${g}S:607 

  guesscoord
  regenerate angles dihedrals
  writepsf "my_4zmj_protomer_${g}${b}_glycans.psf"
  writepdb "my_4zmj_protomer_${g}${b}_glycans_rawloops_tmp1.pdb"
  lappend LOCALFILES my_4zmj_protomer_${g}${b}_glycans.psf
  lappend LOCALFILES my_4zmj_protomer_${g}${b}_glycans_rawloops_tmp1.pdb
}

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

