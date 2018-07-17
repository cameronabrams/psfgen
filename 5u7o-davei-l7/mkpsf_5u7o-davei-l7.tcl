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
set seed 12345
set MPER_EXTEND 0
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "+protomer" } {
     set protomer_only 1
  }
  if { $arg == "-seed" } {
     incr a
     set seed [lindex $argv $a]
  }
  if { $arg == "-mper-extend" } {
     set MPER_EXTEND 1
  }
}

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

set SYSNAME 5u7o-davei-l7

expr srand($seed)

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
set bms529_no_azole [atomselect top "resname 83J and name H3 H14 H6 H5 H4 C23 C26 C24 C21 C20 C15 C13 O06 N05 C10 H1 H2 C07 H11 H12 C04 H9 H10 C01 H7 H8 N02 C12 O03 C14 O09 C16 C18 H13 N08 H23 C19 C17 C22 O11 C27 H16 H17 H18 C25 H15 N28 C29"]
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
$bms529_no_azole set chain "G"
$bms529_no_azole set segname "X0"
$bms529_no_azole set resid 1
$bms529_no_azole set resname AEG
$gp120_glycan set segname "G0g"
$gp41_glycan set segname "B0g"
$fab_pgt122_h_glycan set segname "H0g"
$so4 set segname "S0"

set b529_anlist [list H3 H14 H6 H5 H4 C23 C26 C24 C21 C20 C15 C13 O06 N05 C10 H1 H2 C07 H11 H12 C04 H9 H10 C01 H7 H8 N02 C12 O03 C14 O09 C16 C18 H13 N08 H23 C19 C17 C22 O11 C27 H16 H17 H18 C25 H15 N28 C29]
set aeg_anlist [list H2 H1 H5 H4 H3 C1 C2 C3 C4 C6 C5 C7 O1 N1 C9 H8 H9 C11 H10 H11 C8 H6 H7 C10 H12 H13 N2 C12 O2 C13 O3 C14 C15 H15 N3 H14 C17 C16 C18 O4 C21 H17 H18 H19 C19 H26 N4 C20]
foreach i $b529_anlist j $aeg_anlist {
  [atomselect top "resname AEG and name $i"] set name "x$j"
}
foreach j $aeg_anlist {
  [atomselect top "resname AEG and name x$j"] set name $j
}

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
  $bms529_no_azole set chain "G"
  $bms529_no_azole set segname "X1"
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
  $bms529_no_azole set chain "G"
  $bms529_no_azole set segname "X2"
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
topology $PSFGEN_BASEDIR/charmm/aeg.str
topology $PSFGEN_BASEDIR/charmm/dls1.str
topology $PSFGEN_BASEDIR/charmm/dls2.str
topology $PSFGEN_BASEDIR/charmm/al1p.str

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

  set s10 [atomselect $molid "resname AEG and name H2 H1 H5 H4 H3 C1 C2 C3 C4 C6 C5 C7 O1 N1 C9 H8 H9 C11 H10 H11 C8 H6 H7 C10 H12 H13 N2 C12 O2 C13 O3 C14 C15 H15 N3 H14 C17 C16 C18 O4 C21 H17 H18 H19 C19 H26 N4 C20 C22"]
  $s10 set segname ${x}
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
    if { $MPER_EXTEND == "1" } {
      residue 664 ASP ${b}
      residue 665 LYS ${b}
      residue 666 TRP ${b}
      residue 667 ALA ${b}
      residue 668 SER ${b}
      residue 669 LEU ${b}
      residue 670 TRP ${b}
      residue 671 ASN ${b}
      residue 672 TRP ${b}
    }
  }

  segment ${g}S {
    pdb ${g}-s8.pdb
  }

  segment ${b}S {
    pdb ${b}-s9.pdb
  }

  segment ${x} {
    pdb ${x}-s10.pdb
  }

  segment ${x}1 {
    residue 1 DLS1 ${x}
  }

  segment ${x}2 {
    residue 1 DLS2 ${x}
    residue 2 DLS2 ${x}
    residue 3 DLS2 ${x}
    residue 4 DLS2 ${x}
    residue 5 DLS2 ${x}
    residue 6 DLS2 ${x}
    residue 7 DLS2 ${x}
  }

  segment ${x}T {
    first none
    residue 1 ASP ${x}
    residue 2 LYS ${x}
    residue 3 TRP ${x}
    residue 4 ALA ${x}
    residue 5 SER ${x}
    residue 6 ILE ${x}
    residue 7 TRP ${x}
    residue 8 ASN ${x}
    residue 9 TRP ${x}
    last CT2
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

# bonds in DAVEI
  patch AL1P ${x}:1 ${x}1:1
  patch LL12 ${x}1:1 ${x}2:1
  patch LL22 ${x}2:1 ${x}2:2
  patch LL22 ${x}2:2 ${x}2:3
  patch LL22 ${x}2:3 ${x}2:4
  patch LL22 ${x}2:4 ${x}2:5
  patch LL22 ${x}2:5 ${x}2:6
  patch LL22 ${x}2:6 ${x}2:7
  patch LL2P ${x}2:7 ${x}T:1

  guesscoord
  regenerate angles dihedrals
  writepsf "my_${SYSNAME}_protomer_${g}.psf"
  writepdb "my_${SYSNAME}_protomer_${g}_rawloops_tmp1.pdb"
  lappend LOCALFILES my_${SYSNAME}_protomer_${g}.psf
  lappend LOCALFILES my_${SYSNAME}_protomer_${g}_rawloops_tmp1.pdb
}

resetpsf

foreach g $glist b $blist {
  readpsf my_${SYSNAME}_protomer_${g}.psf pdb my_${SYSNAME}_protomer_${g}_rawloops_tmp1.pdb
}
writepsf "my_${SYSNAME}.psf"
writepdb "my_${SYSNAME}.pdb"
#lappend LOCALFILES my_${SYSNAME}.pdb

mol new my_${SYSNAME}.psf
mol addfile my_${SYSNAME}.pdb
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

$a writepdb "my_${SYSNAME}_fix.pdb"

$a set beta 0

# make a slight manual adjustment to the backbone trajectory for the
# unstructured residues added in gp41 to prevent massive collisions
foreach b $blist {
  set residueList [[atomselect top "chain ${b} and resid 547 to 567 and name CA"] get residue]
  set r1 [lindex $residueList 0]
  set r2 [lindex $residueList end]
  Crot_phi $r1 $r2 ${b} $molid -180
}

foreach x $xlist {
   set sel [atomselect $molid "segname ${x}T"]
   fold_alpha_helix $molid $sel
   $sel delete
}

if { $MPER_EXTEND == "1" } {
   foreach b $blist {
      set sel [atomselect $molid "protein and chain $b and resid 662 to 672"]
      fold_alpha_helix $molid $sel
      $sel delete
   }
}

$a writepdb "my_${SYSNAME}_mcIn.pdb"
lappend LOCALFILES my_${SYSNAME}_mcIn.pdb

set nc 1000
set rcut 3.0
set temperature 2.5
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]
if { 1 } {
foreach g $glist b $blist {
  set residueList [[atomselect ${molid} "chain ${g} and resid 138 to 142 and name CA"] get residue]
  do_loop_mc ${residueList} ${g} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 100 9999]
  set residueList [[atomselect ${molid} "chain ${g} and resid 186 to 189 and name CA"] get residue]
  do_loop_mc ${residueList} ${g} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 100 9999]
  set residueList [[atomselect ${molid} "chain ${g} and resid 399 to 409 and name CA"] get residue]
  do_loop_mc ${residueList} ${g} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 100 9999]
  set residueList [[atomselect ${molid} "chain ${b} and resid 547 to 567 and name CA"] get residue]
  do_loop_mc ${residueList} ${b} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 100 9999]
}
}
set k [expr $k / 10]
set temperature 2.5
set nc 2000
#set bg2 [atomselect top "name CA or (noh and resname AEG DLS1 DLS2)"] 

foreach x $xlist b $blist bir [list 1 2 0] {
   set msel [atomselect ${molid} "(chain $x and resname AEG DLS1 DLS2) or (segname ${x}T)"]
   set ri [[atomselect ${molid} "(chain $x and ((resname AEG and name C24 C26) or (resname DLS1 and name C1 C2 O2 C3 C4 O3 C5 C6 O4 C7 C8 O5 C9 C10 C12 C13 C14) or (resname DLS2 and name N C1 C2 O1 C3 C4 O2 C5))) or (segname ${x}T and resname ASP and resid 1 and name N)"] get index]
   set rj [[atomselect ${molid} "(chain $x and ((resname AEG and name C26) or (resname DLS1 and name N1 C2 O2 C3 C4 O3 C5 C6 O4 C7 C8 O5 C9 C10 N2 C13 C14 C15) or (resname DLS2 and name C1 C2 O1 C3 C4 O2 C5 C6))) or (segname ${x}T and resname ASP and resid 1 and name CA)"] get index]
   set i [[atomselect ${molid} "segname ${x}T and resid 9 and name CA"] get index]
   set j [[atomselect ${molid} "chain $b and resid 641 and name CA"] get index]
   if { $MPER_EXTEND == "1" } {
     set j [[atomselect ${molid} "chain [lindex $blist $bir] and resid 670 and name CA"] get index]
   }
   set fa [[atomselect ${molid} "chain ${x} and resname AEG and name C1"] get index]
   do_flex_mc $molid $msel $ri $rj $fa $k $i $j $bg $rcut $nc $temperature [irand_dom 100 9999]
}

$a writepdb "my_${SYSNAME}_mcOut.pdb"

foreach f $LOCALFILES {
  exec rm $f
}

quit

