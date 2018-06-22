# VMD/psfgen script for generating psf/pdb pair for PDB 5u7o
# trimeric HIV-1 gp140 with BMS529 bound
#
# cameron f abrams (c) 2018
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

# check for any arguments
set protomer_only 0
#set protomer_only 1
foreach arg $argv {
  if { $arg == "+protomer" } {
     set protomer_only 1
  }
}

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

mol new 5u7o.pdb
set a [atomselect top all]
set g [atomselect top "chain G and name CA"]
# align the asymmetric unit onto 4TVP, so we can use 4TVP's BIOMT
# transformations to construct the trimer (since the 5U7O entry did not
# contain any).
mol new 4tvp.pdb
set g0 [atomselect top "chain G and name CA and resid 31 to 60 62 to 137 151 to 185 190 to 505"]
$a move [measure fit $g $g0]
mol delete top

# select separate parts
set gp120 [atomselect top "protein and chain G"]
set gp41 [atomselect top "protein and chain B"]
set fab_35o22_h [atomselect top "protein and chain D"]
set fab_35o22_l [atomselect top "protein and chain E"]
set fab_pgt122_h [atomselect top "protein and chain H"]
set fab_pgt122_l [atomselect top "protein and chain L"]
set bms529 [atomselect top "resname 83J"]
set gp120_glycan [atomselect top "chain G and resname NAG BMA MAN"]
set gp41_glycan [atomselect top  "chain B and resname NAG BMA MAN"]
set fab_pgt122_h_glycan [atomselect top "chain H and resname NAG"]
set so4 [atomselect top "resname SO4"]

# update segnames and save coordinates
$gp120 set segname "G0"
$gp41 set segname "B0"
$fab_35o22_h set segname "D0"
$fab_35o22_l set segname "E0"
$fab_pgt122_h set segname "H0"
$fab_pgt122_l set segname "L0"
$bms529 set chain "X"
$bms529 set segname "X0"
$gp120_glycan set segname "G0g"
$gp41_glycan set segname "B0g"
$fab_pgt122_h_glycan set segname "H0g"
$so4 set segname "S0"

$a writepdb "unit0.pdb"

if { $protomer_only == "0"} {

  set xx [list -0.500000 -0.866025  0.000000  -515.56   ]
  set yy [list  0.866025 -0.500000  0.000000     0      ]
  set zz [list  0.000000  0.000000  1.000000     0.00000]
  set tt [list  0         0         0            1      ]
  set tmat3 [list $xx $yy $zz $tt]

  $a move $tmat3
  $gp120 set segname "G1"
  $gp41 set segname "B1"
  $fab_35o22_h set segname "D1"
  $fab_35o22_l set segname "E1"
  $fab_pgt122_h set segname "H1"
  $fab_pgt122_l set segname "L1"
  $bms529 set chain "Y"
  $bms529 set segname "X1"
  $gp120_glycan set segname "G1g"
  $gp41_glycan set segname "B1g"
  $fab_pgt122_h_glycan set segname "H1g"
  $so4 set segname "S1"

  $a writepdb "unit1.pdb"; lappend LOCALFILES unit1.pdb

  $a move $tmat3
  $gp120 set segname "G2"
  $gp41 set segname "B2"
  $fab_35o22_h set segname "D2"
  $fab_35o22_l set segname "E2"
  $fab_pgt122_h set segname "H2"
  $fab_pgt122_l set segname "L2"
  $bms529 set chain "Z"
  $bms529 set segname "X2"
  $gp120_glycan set segname "G2g"
  $gp41_glycan set segname "B2g"
  $fab_pgt122_h_glycan set segname "H2g"
  $so4 set segname "S2"

  $a writepdb "unit2.pdb"; lappend LOCALFILES unit2.pdb
}

mol delete top

package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
#topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology $PSFGEN_BASEDIR/charmm/bms529.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias atom BGNA N2 N
pdbalias atom BGNA C7 C
pdbalias atom BGNA O7 O
pdbalias atom BGNA C8 CT

# chain mappings
# A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
#   x   o @   ! $       *                       #      unit 0
#     x     o       @ !   $ *                     #    unit 1
# x                             o @ ! $ *           #  unit 2
set ulist [list 0 1 2]
set glist [list G K R]
set blist [list B C A]
set dlist [list D F P]
set elist [list E J Q]
set hlist [list H M R]
set llist [list L N T]
set xlist [list X Y Z]
if { $protomer_only == "1" } {
  set ulist [list 0]
  set glist [list G]
  set blist [list B]
  set dlist [list D]
  set elist [list E]
  set hlist [list H]
  set llist [list L]
  set xlist [list X]
}

foreach u $ulist g $glist b $blist d $dlist e $elist h $hlist l $llist x $xlist {
 
  mol new "unit${u}.pdb"
  set molid [molinfo top get id]

  set s1 [atomselect $molid "protein and chain G and resid 31 to 60"]
  $s1 set chain ${g}
  $s1 writepdb "${g}-s1.pdb"; lappend LOCALFILES ${g}-s1.pdb
  set s1n_pos [cacoIn_nOut  60 ${g} $molid]
  set s2 [atomselect $molid "protein and chain G and resid 62 to 137"]
  $s2 set chain ${g}
  $s2 writepdb "${g}-s2.pdb"; lappend LOCALFILES ${g}-s2.pdb
  set s2n_pos [cacoIn_nOut 137 ${g} $molid]
  set s3 [atomselect $molid "protein and chain G and resid 151 to 185"]
  $s3 set chain ${g}
  $s3 writepdb "${g}-s3.pdb"; lappend LOCALFILES ${g}-s3.pdb
  set s3n_pos [cacoIn_nOut 185 ${g} $molid]
  set s4 [atomselect $molid "protein and chain G and resid 190 to 398"]
  $s4 set chain ${g}
  $s4 writepdb "${g}-s4.pdb"; lappend LOCALFILES ${g}-s4.pdb
  set s4n_pos [cacoIn_nOut 398 ${g} $molid]
  set s5 [atomselect $molid "protein and chain G and resid 411 to 505"]
  $s5 set chain ${g}
  $s5 writepdb "${g}-s5.pdb"; lappend LOCALFILES ${g}-s5.pdb

  set s6 [atomselect $molid "protein and chain B and resid 518 to 547"]
  $s6 set chain $b
  $s6 writepdb "${b}-s6.pdb"; lappend LOCALFILES ${b}-s6.pdb
  set s6n_pos [cacoIn_nOut 547 $b $molid]
  set s7 [atomselect $molid "protein and chain B and resid 568 to 663"]
  $s7 set chain $b
  $s7 writepdb "${b}-s7.pdb"; lappend LOCALFILES ${b}-s7.pdb

  [atomselect $molid "resname NAG"] set resname BGNA
  [atomselect $molid "resname MAN"] set resname AMAN
  [atomselect $molid "resname BMA"] set resname BMAN
  set s8 [atomselect $molid "chain G and resname BGNA AMAN BMAN"]
  $s8 set chain $g
  $s8 set segname ${g}S
  $s8 writepdb "${g}-s8.pdb"; lappend LOCALFILES ${g}-s8.pdb
  set s9 [atomselect $molid "chain B and resname BGNA AMAN BMAN"]
  $s9 set chain $b
  $s9 set segname ${b}S
  $s9 writepdb "${b}-s9.pdb"; lappend LOCALFILES ${b}-s9.pdb

  set s10 [atomselect $molid "resname 83J"]
  $s10 set resname B529
  $s10 set segname $x
  $s10 set chain $x
  $s10 writepdb "${x}-s10.pdb"; lappend LOCALFILES ${x}-s10.pdb
  
  mol delete $molid

  resetpsf

  segment ${g} {
    pdb ${g}-s1.pdb
    residue 61 TYR ${g}
    pdb ${g}-s2.pdb
    residue 138 ILE ${g}
    residue 139 THR ${g}
    residue 140 ASP ${g}
    residue 141 ASP ${g}
    residue 142 MET ${g}
    pdb ${g}-s3.pdb
    residue 186 GLU ${g}
    residue 187 ASN ${g}
    residue 188 GLN ${g}
    residue 189A GLY ${g}
    residue 189B ASN ${g}
    residue 189C ARG ${g}
    residue 189D SER ${g}
    residue 189E ASN ${g}
    residue 189F ASN ${g}
    residue 189G SER ${g}
    residue 189H ASN ${g}
    residue 189I LYS ${g}
    pdb ${g}-s4.pdb
    residue 399 THR ${g}
    residue 400 SER ${g}
    residue 401 VAL ${g}
    residue 402 GLN ${g}
    residue 403 GLY ${g}
    residue 404 SER ${g}
    residue 405 ASN ${g}
    residue 406 SER ${g}
    residue 407 THR ${g}
    residue 408 GLY ${g}
    residue 409 SER ${g}
    pdb ${g}-s5.pdb
  }

  segment ${b} {
    residue 512 ALA ${b}
    residue 513 VAL ${b}
    residue 514 GLY ${b}
    residue 515 ILE ${b}
    residue 516 GLY ${b}
    residue 517 ALA ${b}
    pdb ${b}-s6.pdb
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
    pdb ${b}-s7.pdb 
  }

  segment ${g}S {
    pdb ${g}-s8.pdb
  }

  segment ${b}S {
    pdb ${b}-s9.pdb
  }

  segment $x {
    pdb ${x}-s10.pdb
  }

  coordpdb ${g}-s1.pdb ${g}
  coordpdb ${g}-s2.pdb ${g}
  coordpdb ${g}-s3.pdb ${g}
  coordpdb ${g}-s4.pdb ${g}
  coordpdb ${g}-s5.pdb ${g}
  coordpdb ${b}-s6.pdb ${b}
  coordpdb ${b}-s7.pdb ${b}
  coordpdb ${g}-s8.pdb ${g}S
  coordpdb ${b}-s9.pdb ${b}S
  coordpdb ${x}-s10.pdb ${x}

  coord ${g} 61  N $s1n_pos
  coord ${g} 138 N $s2n_pos
  coord ${g} 186 N $s3n_pos
  coord ${g} 399 N $s4n_pos
  coord ${b} 548 N $s6n_pos

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

# Glycan LINK's on chains G (55) and B (3) (58 total)
# NGLB patches: 18 on chain G, 3 on chain B (21 total so far)
#   grep ^LINK 5u7o.pdb | grep ND2 | grep -w G | awk '{printf("  patch NGLB ${g}:%i ${g}S:%i\n",$5,$9);}'
  patch NGLB ${g}:88 ${g}S:604
  patch NGLB ${g}:133 ${g}S:611
  patch NGLB ${g}:137 ${g}S:612
  patch NGLB ${g}:156 ${g}S:616
  patch NGLB ${g}:160 ${g}S:621
  patch NGLB ${g}:197 ${g}S:623
  patch NGLB ${g}:234 ${g}S:625
  patch NGLB ${g}:262 ${g}S:627
  patch NGLB ${g}:276 ${g}S:633
  patch NGLB ${g}:295 ${g}S:634
  patch NGLB ${g}:301 ${g}S:636
  patch NGLB ${g}:332 ${g}S:638
  patch NGLB ${g}:339 ${g}S:658
  patch NGLB ${g}:355 ${g}S:648
  patch NGLB ${g}:363 ${g}S:649
  patch NGLB ${g}:386 ${g}S:651
  patch NGLB ${g}:392 ${g}S:653
  patch NGLB ${g}:448 ${g}S:655
#   grep ^LINK 5u7o.pdb | grep ND2 | grep -w B | awk '{printf("  patch NGLB ${b}:%i ${b}S:%i\n",$5,$9);}'
  patch NGLB ${b}:611 ${b}S:702
  patch NGLB ${b}:618 ${b}S:703
  patch NGLB ${b}:637 ${b}S:704

# 12aa patches: 4 on chain G, 0 on chain B  (25 total so far)
#   grep ^LINK 5u7o.pdb | grep -w O2 | grep -w G | awk '{printf("  patch 12aa ${g}S:%i ${g}S:%i\n",$5,$9);}'
  patch 12aa ${g}S:631 ${g}S:632
  patch 12aa ${g}S:643 ${g}S:647
  patch 12aa ${g}S:644 ${g}S:645
  patch 12aa ${g}S:645 ${g}S:646

# 13ab patches: 7 on chain G, 0 on chain B (32 total so far)
#   grep ^LINK 5u7o.pdb | grep -w O3 | awk '{printf("  patch 13ab ${g}S:%i ${g}S:%i\n",$5,$9);}'
  patch 13ab ${g}S:606 ${g}S:607
  patch 13ab ${g}S:608 ${g}S:610
  patch 13ab ${g}S:614 ${g}S:615
  patch 13ab ${g}S:618 ${g}S:619
  patch 13ab ${g}S:629 ${g}S:631
  patch 13ab ${g}S:640 ${g}S:644
  patch 13ab ${g}S:641 ${g}S:643

# 14bb patches: 20 on chain G, 0 on chain B (52 total so far)
#   grep ^LINK 5u7o.pdb | grep -w O4 | awk '{printf("  patch 14bb ${g}S:%i ${g}S:%i\n",$5,$9);}'
  patch 14bb ${g}S:604 ${g}S:605
  patch 14bb ${g}S:605 ${g}S:606
  patch 14bb ${g}S:612 ${g}S:613
  patch 14bb ${g}S:613 ${g}S:614
  patch 14bb ${g}S:616 ${g}S:617
  patch 14bb ${g}S:617 ${g}S:618
  patch 14bb ${g}S:621 ${g}S:622
  patch 14bb ${g}S:623 ${g}S:624
  patch 14bb ${g}S:625 ${g}S:626
  patch 14bb ${g}S:627 ${g}S:628
  patch 14bb ${g}S:628 ${g}S:629
  patch 14bb ${g}S:634 ${g}S:635
  patch 14bb ${g}S:636 ${g}S:637
  patch 14bb ${g}S:638 ${g}S:639
  patch 14bb ${g}S:639 ${g}S:640
  patch 14bb ${g}S:649 ${g}S:650
  patch 14bb ${g}S:651 ${g}S:652
  patch 14bb ${g}S:653 ${g}S:654
  patch 14bb ${g}S:655 ${g}S:656
  patch 14bb ${g}S:656 ${g}S:657

# 16ab patches: 6 on chain G, 0 on chain B (58 total)
#   grep ^LINK 5u7o.pdb | grep -w O6 | awk '{printf("  patch 16ab ${g}S:%i ${g}S:%i\n",$5,$9);}'
  patch 16ab ${g}S:606 ${g}S:608
  patch 16ab ${g}S:608 ${g}S:609
  patch 16ab ${g}S:618 ${g}S:620
  patch 16ab ${g}S:629 ${g}S:630
  patch 16ab ${g}S:640 ${g}S:641
  patch 16ab ${g}S:641 ${g}S:642

  guesscoord
  regenerate angles dihedrals
  writepsf "my_5u7o_protomer_${g}.psf"
  writepdb "my_5u7o_protomer_${g}_rawloops_tmp1.pdb"
  lappend LOCALFILES my_5u7o_protomer_${g}.psf
  lappend LOCALFILES my_5u7o_protomer_${g}_rawloops_tmp1.pdb
}

resetpsf

foreach g $glist b $blist {
  readpsf my_5u7o_protomer_${g}.psf pdb my_5u7o_protomer_${g}_rawloops_tmp1.pdb
}
writepsf "my_5u7o.psf"
writepdb "my_5u7o.pdb"
#lappend LOCALFILES my_5u7o.pdb

mol new my_5u7o.psf
mol addfile my_5u7o.pdb
set molid [molinfo top get id]

set a [atomselect top all]
$a set beta 1 ;# pretty much everything is fixed except...

foreach g $glist b $blist {
  [atomselect top "chain $g and resid 61"] set beta 0
  [atomselect top "chain $g and resid 138 to 142"] set beta 0
  [atomselect top "chain $g and resid 186 to 189"] set beta 0
  [atomselect top "chain $g and resid 399 to 409"] set beta 0
  [atomselect top "chain $b and resid 512 to 517"] set beta 0
  [atomselect top "chain $b and resid 548 to 567"] set beta 0
} 

$a writepdb "my_5u7o_fix.pdb"

$a set beta 0

# make a slight manual adjustment to the backbone trajectory for the
# unstructured residues added in gp41 to prevent massive collisions
foreach b $blist {
  set residueList [[atomselect top "chain ${b} and resid 547 to 567 and name CA"] get residue]
  set r1 [lindex $residueList 0]
  set r2 [lindex $residueList end]
  Crot_phi $r1 $r2 ${b} $molid -180
}

$a writepdb "my_5u7o_mcIn.pdb"
#lappend LOCALFILES my_5u7o_mcIn.pdb

set nc 1000
set rcut 3.0
set temperature 2.5
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]

foreach g $glist b $blist {
  set residueList [[atomselect ${molid} "chain ${g} and resid 186 to 189 and name CA"] get residue]
  do_loop_mc ${residueList} ${g} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 100 9999]
  set residueList [[atomselect ${molid} "chain ${g} and resid 399 to 409 and name CA"] get residue]
  do_loop_mc ${residueList} ${g} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 100 9999]
  set residueList [[atomselect ${molid} "chain ${b} and resid 547 to 567 and name CA"] get residue]
  do_loop_mc ${residueList} ${b} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 100 9999]
}

$a writepdb "my_5u7o_mcOut.pdb"

foreach f $LOCALFILES {
  exec rm $f
}

quit

