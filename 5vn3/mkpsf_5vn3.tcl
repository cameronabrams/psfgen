# VMD/psfgen script for generating psf/pdb pair for PDB 5vn3
# trimeric HIV-1 gp140 in "open" state
# 
# note: 17b antibodies (6 chains) and sCD4's (3 chains) are 
# not included as yet.
#
# cameron f abrams (c) 2017-18
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
set WARHEAD_SEQUENCE [list ASP LYS TRP ALA SER ILE TRP ASN TRP ]
set MPER_665_to_682 [list LYS TRP ALA SER LEU TRP ASN TRP PHE ASP ILE SER ASN TRP LEU TRP TYR ILE]
set MPER_EXTEND 0
set seed 12345
set attractor_radius 60
set SKIP_LOOPMC 0
set LOG_DCD 0
set logid -1
set DOCK_BNM 0
set DAVEI 0
set LOOP_MC_NC  1000
set LOOP_MC_RCUT  3.0
set LOOP_MC_TEMPERATURE 2.5
set LOOP_MC_K  1.0
set LOOP_MC_R0 2.0
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-seed" } {
     incr a
     set seed [lindex $argv $a]
  }
  if { $arg == "-mper-extend" } {
     set MPER_EXTEND 1
  }
  if { $arg == "-dock-bnm" } {
     set DOCK_BNM 1
  }
  if { $arg == "-skip-loopmc" } {
     set SKIP_LOOPMC 1
  }
  if { $arg == "-log-dcd" } {
     set LOG_DCD 1
     incr a
     set log_dcd_file [lindex $argv $a]
  }
  if { $arg == "-loop-mc-numcycles" || $arg == "-lmcnc" } {
     incr a
     set LOOP_MC_NC [lindex $argv $a]
  }
  if { $arg == "-loop-mc-cutoff" || $arg == "-lmcrc" } {
     incr a
     set LOOP_MC_RCUT [lindex $argv $a]
  }
  if { $arg == "-loop-mc-temperature" || $arg == "-lmct" } {
     incr a
     set LOOP_MC_TEMPERATURE [lindex $argv $a]
  }
  if { $arg == "-loop-mc-k" || $arg == "-lmck" } {
     incr a
     set LOOP_MC_K [lindex $argv $a]
  }
  if { $arg == "-loop-mc-r0" || $arg == "-lmcr0" } {
     incr a
     set LOOP_MC_R0 [lindex $argv $a]
  }
  if { $arg == "-davei" } {
     incr a
     set DAVEI [lindex $argv $a]
  }
  if { $arg == "-attractor_radius" } {
     incr a
     set attractor_radius [lindex $argv $a]
  }
}

if { [expr $DAVEI > 0] == "1" } {
   set DOCK_BNM 1
}

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

mol new 5vn3.pdb
set molid [molinfo top get id]

set glist { G I J }
set blist { A B D }
set xlist { X Y Z }

set gsegs {
  {{32 131} {193 309} {325 404} {412 505}}
  {{32 131} {193 309} {325 404} {412 505}}
  {{32 131} {193 309} {325 404} {412 505}} }
set bsegs {
  {{515 548} {563 664}}
  {{515 548} {563 664}}
  {{515 548} {563 664}} }

set ggaps {
  {{  132 134 { GLY GLY GLY } }
   {  311 324 { GLY PRO GLY ARG ALA PHE TYR ALA THR GLY ASP ILE ILE GLY } }
   {  407 411 { PRO THR GLY GLY GLU } }}
  {{  132 134 { GLY GLY GLY } }
   {  311 324 { GLY PRO GLY ARG ALA PHE TYR ALA THR GLY ASP ILE ILE GLY } }
   {  407 411 { PRO THR GLY GLY GLU } }}
  {{  132 134 { GLY GLY GLY } }
   {  311 324 { GLY PRO GLY ARG ALA PHE TYR ALA THR GLY ASP ILE ILE GLY } }
   {  407 411 { PRO THR GLY GLY GLU } }}}
set bgaps {
  {{  549 562 { VAL GLN GLN GLN ASN ASN LEU LEU ARG ALA PRO GLU ALA GLN }}}
  {{  549 562 { VAL GLN GLN GLN ASN ASN LEU LEU ARG ALA PRO GLU ALA GLN }}}
  {{  549 562 { VAL GLN GLN GLN ASN ASN LEU LEU ARG ALA PRO GLU ALA GLN }}}}
set ssbonds {
 { G  54 G  74 } { G 119 G 205 } { G 126 G 196 } { G 218 G 247 } 
 { G 228 G 239 } { G 296 G 331 } { G 378 G 445 } { G 385 G 418 } 
 { A 598 A 604 } { A 605 G 501 }
 { I  54 I  74 } { I 119 I 205 } { I 126 I 196 } { I 218 I 247 } 
 { I 228 I 239 } { I 296 I 331 } { I 378 I 445 } { I 385 I 418 } 
 { B 598 B 604 } { B 605 I 501 }
 { J  54 J  74 } { J 119 J 205 } { J 126 J 196 } { J 218 J 247 } 
 { J 228 J 239 } { J 296 J 331 } { J 378 J 445 } { J 385 J 418 } 
 { D 598 D 604 } { D 605 J 501 }
}

set mln {}
foreach g $glist b $blist {
  lappend mln [list $g 132 131]
  lappend mln [list $g 311 309]
  lappend mln [list $g 407 404]
  lappend mln [list $b 549 548]
} 
set ns {}
foreach m $mln {
  lappend ns [cacoIn_nOut [lindex $m 2] [lindex $m 0] $molid]
}

for { set m 0 } { $m < 3 } { incr m } {
   set g [lindex $glist $m]
   set gs [lindex $gsegs $m]
   set ngs [llength $gs]
   for { set s 0 } { $s < $ngs } { incr s } {
     set ts [lindex $gs $s]
     set nr [lindex $ts 0]
     set cr [lindex $ts 1]
     [atomselect $molid "protein and chain $g and resid $nr to $cr"] writepdb "${g}_${nr}_to_${cr}.pdb"
     lappend LOCALFILES ${g}_${nr}_to_${cr}.pdb
   }
   set b [lindex $blist $m]
   set bs [lindex $bsegs $m]
   set nbs [llength $bs]
   for { set s 0 } { $s < $nbs } { incr s } {
     set ts [lindex $bs $s]
     set nr [lindex $ts 0]
     set cr [lindex $ts 1]
     [atomselect $molid "protein and chain $b and resid $nr to $cr"] writepdb "${b}_${nr}_to_${cr}.pdb"
     lappend LOCALFILES ${b}_${nr}_to_${cr}.pdb
   }   

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

if { $DOCK_BNM == "1" } {
  mol new 5f4p.pdb
  set refid [molinfo top get id]
  mol top $molid

  set bnm [atomselect $refid "chain A and resname 5VG"]
  $bnm set resid 1
  $bnm set resname BNM3
  set origcoor [$bnm get {x y z}]

  set pocketresidlist [[atomselect $refid "name CA and (same residue as exwithin 5.0 of resname 5VG)"] get resid]
  set refpocket [atomselect $refid "chain A and resid $pocketresidlist and name CA C N"]

  for { set m 0 } { $m < 3 } { incr m } {
    set x [lindex $xlist $m]
    set g [lindex $glist $m]
    $bnm set chain $x
    # align, move
    set targpocket [atomselect $molid "chain $g and resid $pocketresidlist and name CA C N"]
    if { [$refpocket num] != [$targpocket num] } {
      puts "ERROR: bnm-pockets not congruent"
      exit
    }
    $bnm move [measure fit $refpocket $targpocket]
    $bnm writepdb "${x}.pdb"
    $bnm set {x y z} $origcoor
  }

  mol delete $refid
}

mol delete $molid

package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology ${PSFGEN_BASEDIR}/charmm/bnm.str
if { [expr $DAVEI > 0] == "1" } {
  topology ${PSFGEN_BASEDIR}/charmm/dls1.str
  topology ${PSFGEN_BASEDIR}/charmm/dls2.str
  topology ${PSFGEN_BASEDIR}/charmm/al1p_PEG.str
}

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias atom BGNA N2 N
pdbalias atom BGNA C7 C
pdbalias atom BGNA O7 O
pdbalias atom BGNA C8 CT
foreach nbad [exec grep 5VG 5f4p.pdb | grep -w "5VG A" | grep HETATM | cut -b 13-16] ngood [exec grep ATOM ${PSFGEN_BASEDIR}/charmm/bnm.str | cut -b 6-9 | grep -v ^H] {
  pdbalias atom BNM3 $nbad $ngood
}

for { set m 0 } { $m < 3 } { incr m } {
  set g [lindex $glist $m]
  set b [lindex $blist $m]
  set x [lindex $xlist $m]

  set gs [lindex $gsegs $m]
  set ngs [llength $gs]
  set gg [lindex $ggaps $m]
  set ngg [llength $gg]

  set bs [lindex $bsegs $m]
  set nbs [llength $bs]
  set bg [lindex $bgaps $m]
  set nbg [llength $bg]

  segment $g {
    for { set s 0 } { $s < $ngs } { incr s } {
       set ts [lindex $gs $s]
       set nr [lindex $ts 0]
       set cr [lindex $ts 1]
       pdb ${g}_${nr}_to_${cr}.pdb
       if { $s < $ngg } {
          set tg [lindex $gg $s]
          set nrr [lindex $tg 0]
          set crr [lindex $tg 1]
          set seq [lindex $tg 2]
          foreach r $seq {
            residue $nrr $r $g
            set nrr [my_increment $nrr]
          }
       }
    }
  }

  segment ${g}S {
    pdb ${g}S.pdb
  }

  segment $b {
    for { set s 0 } { $s < $nbs } { incr s } {
       set ts [lindex $bs $s]
       set nr [lindex $ts 0]
       set cr [lindex $ts 1]
       pdb ${b}_${nr}_to_${cr}.pdb
       if { $s < $nbg } {
          set tg [lindex $bg $s]
          set nrr [lindex $tg 0]
          set crr [lindex $tg 1]
          set seq [lindex $tg 2]
          foreach r $seq {
            residue $nrr $r $b
            set nrr [my_increment $nrr]
          }
       }
    }
    if { $MPER_EXTEND == "1" } {
      set lr 0
      for { set r 665 } { $r < 683 } { incr r } {
        residue $r [lindex $MPER_665_to_682 $lr] $b
        incr lr
      }
    }
  }
  
  segment ${b}S {
    pdb ${b}S.pdb
  }

  if { $DOCK_BNM == "1" } {
#    set x [lindex $xlist $m]
    segment ${x} {
      pdb ${x}.pdb
    }
    if { [expr $DAVEI > 0] == "1" } {
       # DLS1 linker (PEG/azide)
       segment ${x}1 {
         residue 1 DLS1 ${x}
       }
       # DLS2 aminodiethoxyacetate units
       segment ${x}2 {
         for { set ll 0 } { $ll < $DAVEI } { incr ll } {
            set lll [ expr $ll + 1 ]
            residue $lll DLS2 ${x}
         }
       }
       # Warhead peptide (e.g., Trp3) with no N-terminal patch and a neutral C-terminal patch
       segment ${x}T {
         first none
         for { set ll 0 } { $ll <  [llength $WARHEAD_SEQUENCE] } { incr ll } {
           set lll [expr $ll + 1]
           residue $lll [lindex $WARHEAD_SEQUENCE $ll] ${x}
         }
         last CT2
       }
    } 
  }
}

for { set m 0 } { $m < 3 } { incr m } {
   set g [lindex $glist $m]
   set gs [lindex $gsegs $m]
   set ngs [llength $gs]
   for { set s 0 } { $s < $ngs } { incr s } {
     set ts [lindex $gs $s]
     set nr [lindex $ts 0]
     set cr [lindex $ts 1]
     coordpdb ${g}_${nr}_to_${cr}.pdb $g
   }
   set b [lindex $blist $m]
   set bs [lindex $bsegs $m]
   set nbs [llength $bs]
   for { set s 0 } { $s < $nbs } { incr s } {
     set ts [lindex $bs $s]
     set nr [lindex $ts 0]
     set cr [lindex $ts 1]
     coordpdb ${b}_${nr}_to_${cr}.pdb $b
   }
   if { $DOCK_BNM == "1" } {
     set x [lindex $xlist $m]
     coordpdb ${x}.pdb ${x}
   }
}
# coordpdb sugars
for { set m 0 } { $m < 3 } { incr m } {
  set g [lindex $glist $m]
  set b [lindex $blist $m]
  coordpdb ${g}S.pdb ${g}S
  coordpdb ${b}S.pdb ${b}S
}

# set N positions
foreach m $mln n $ns {
  puts "coord [lindex $m 0] [lindex $m 1] N $n"
  coord [lindex $m 0] [lindex $m 1] N $n
}

# disulfides
foreach ss $ssbonds {
   set c1 [lindex $ss 0]
   set r1 [lindex $ss 1]
   set c2 [lindex $ss 2]
   set r2 [lindex $ss 3]
   patch DISU ${c1}:${r1} ${c2}:${r2}
}

for { set m 0 } { $m < 3 } { incr m } {
  set g [lindex $glist $m]
  set b [lindex $blist $m]
  set x [lindex $xlist $m]

  patch NGLB ${b}:611 ${b}S:701
  patch NGLB ${b}:616 ${b}S:702
  patch NGLB ${b}:625 ${b}S:703
  patch NGLB ${b}:637 ${b}S:704
  patch NGLB ${g}:88 ${g}S:601
  patch NGLB ${g}:234 ${g}S:602
  patch 14bb ${g}S:602 ${g}S:603
  patch NGLB ${g}:241 ${g}S:604
  patch 14bb ${g}S:604 ${g}S:605
  patch 14bb ${g}S:605 ${g}S:606
  patch NGLB ${g}:262 ${g}S:607
  patch 14bb ${g}S:607 ${g}S:608
  patch 14bb ${g}S:608 ${g}S:609
  patch 13ab ${g}S:609 ${g}S:610
  patch 16ab ${g}S:609 ${g}S:613
  patch 12aa ${g}S:610 ${g}S:611
  patch 12aa ${g}S:611 ${g}S:612
  patch NGLB ${g}:276 ${g}S:614
  patch 14bb ${g}S:614 ${g}S:615
  patch 14bb ${g}S:615 ${g}S:616
  patch 16bb ${g}S:616 ${g}S:617
  patch NGLB ${g}:295 ${g}S:618
  patch 14bb ${g}S:618 ${g}S:619
  patch NGLB ${g}:332 ${g}S:620
  patch 14bb ${g}S:620 ${g}S:621
  patch NGLB ${g}:339 ${g}S:622
  patch NGLB ${g}:355 ${g}S:623
  patch NGLB ${g}:362 ${g}S:624
  patch 14bb ${g}S:624 ${g}S:625
  patch 14bb ${g}S:625 ${g}S:626
  patch 13ab ${g}S:626 ${g}S:627
  patch NGLB ${g}:386 ${g}S:628
  patch 14bb ${g}S:628 ${g}S:629
  patch 14bb ${g}S:629 ${g}S:630
  patch NGLB ${g}:392 ${g}S:637
  patch 14bb ${g}S:637 ${g}S:638
  patch 14bb ${g}S:638 ${g}S:639
  patch NGLB ${g}:397 ${g}S:631
  patch NGLB ${g}:413 ${g}S:632
  patch 14bb ${g}S:632 ${g}S:633
  patch 14bb ${g}S:633 ${g}S:634
  patch NGLB ${g}:448 ${g}S:635
  patch 14bb ${g}S:635 ${g}S:636

  # bonds in DAVEI
  if { [expr $DAVEI > 0] == "1" } {
    patch AL1P ${x}:1 ${x}1:1
    patch LL12 ${x}1:1 ${x}2:1
    for { set nll 1 } { $nll < $DAVEI } { incr nll } {
       patch LL22 ${x}2:${nll} ${x}2:[expr $nll+1]
    }
    patch LL2P ${x}2:${DAVEI} ${x}T:1
  }
  
}

guesscoord
regenerate angles dihedrals

writepsf "my_5vn3.psf"
writepdb "unrelaxed.pdb"

resetpsf

lappend LOCALFILES unrelaxed.pdb

mol delete top
mol new my_5vn3.psf
mol addfile unrelaxed.pdb
set molid [molinfo top get id]

# initialize a logging molecule
if { $LOG_DCD == "1" } {
   # set logid [log_initialize my_${SYSNAME}.psf my_${SYSNAME}.pdb]
   mol new my_5vn3.psf
   mol addfile unrelaxed.pdb
   set logid [molinfo top get id]
   puts "Logging molecule number $logid to be saved to $log_dcd_file"
   mol top $molid
}
# fold MPER extensions into alpha helices
if { $MPER_EXTEND == "1" } {
   for { set m 0 } { $m < 3 } { incr m } {
     set b [lindex $blist $m]
     set Cterm 682
     set sel [atomselect $molid "protein and chain $b and resid 663 to $Cterm"]
     puts "Fold-alpha: chain $b MPER extension"
     fold_alpha_helix $molid $sel extrabonds-${b}.inp
     log_addframe $molid $logid
     $sel delete
   }
}
# fold DAVEI warhead peptides into alpha helices
if { [expr $DAVEI > 0] == "1" } {
  foreach x $xlist {
    puts "Fold-alpha: seg ${x}T"
    set sel [atomselect ${molid} "segname ${x}T"]
    puts "number of atoms in segname ${x}T: [$sel num]"
    fold_alpha_helix $molid $sel extrabonds-${x}T.inp
    lappend LOCALFILES extrabonds-${x}T.inp
    $sel delete
    log_addframe $molid $logid
  }
  exec cat extrabonds-XT.inp extrabonds-YT.inp extrabonds-ZT.inp > extrabonds-TRP3.inp
}

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
log_addframe $molid $logid
if { $logid != "-1" } {
  set la [atomselect $logid all]
  for {set li 0} {$li < [molinfo $logid get numframes]} {incr li} {
    $la frame $li
    $la moveby [vecscale -1 $or]
    $la move [transaxis z [expr -1 * $p] rad]
    $la move [transaxis y [expr -1 * $t] rad]
  }
  $la delete
}

if { $SKIP_LOOPMC == "0" } {
 set nc $LOOP_MC_NC
 set rcut $LOOP_MC_RCUT
 set temperature $LOOP_MC_TEMPERATURE
 set k $LOOP_MC_K
 set r0 $LOOP_MC_R0
 set background [atomselect ${molid} "noh"]

 for { set m 0 } { $m < 3 } { incr m } {
    set g [lindex $glist $m]
    set gg [lindex $ggaps $m]
    set ngg [llength $gg]
    for { set s 0 } { $s < $ngg } { incr s } {
       set tg [lindex $gg $s]
       set nrr [lindex $tg 0]
       set crr [lindex $tg 1]
       set residueList [[atomselect ${molid} "chain $g and resid $nrr to $crr and name CA"] get residue]
       puts "loop mc on chain $g res $nrr to $crr : residueList $residueList"
#       puts "STUB"
       do_loop_mc ${residueList} $g ${molid} ${k} ${r0} ${background} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] $logid
    }
    set b [lindex $blist $m]
    set bg [lindex $bgaps $m]
    set nbg [llength $bg]
    for { set s 0 } { $s < $nbg } { incr s } {
       set tg [lindex $bg $s]
       set nrr [lindex $tg 0]
       set crr [lindex $tg 1]
       set residueList [[atomselect ${molid} "chain $b and resid $nrr to $crr and name CA"] get residue]
       puts "loop mc on chain $b res $nrr to $crr : residueList $residueList"
#       puts "STUB"
       do_loop_mc ${residueList} $b ${molid} ${k} ${r0} ${background} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] $logid
    }
  }
#  puts "HERE"
  if { [ expr $DAVEI > 0 ] == "1" } {
    set k [ expr $k/10 ]
    set out_targ_points {}
    # flexMC to position the ends of the DAVEI's (experimental)
    foreach x $xlist b $blist bir [list 1 2 0] g $glist {
     # 1. identify the flexible chain
     set msel [atomselect ${molid} "chain $x"]
     # 2. identify atoms on the flexible chain that can be "left-hand-sides" of rotatable bonds
     set ri [[atomselect ${molid} "(chain $x and ((resname BNM3 and name C18 C11 N1) or (resname DLS1 and name C3 C4 O3 C5 C6 O4 C7 C8 O5 C9 C10 C12 C13 C14) or (resname DLS2 and name N C1 C2 O1 C3 C4 O2 C5))) or (segname ${x}T and resname ASP and resid 1 and name N)"] get index]
     # 3. identify atoms on the flexible chain that can be "right-hand-sides" of rotatable bonds
     set rj [[atomselect ${molid} "(chain $x and ((resname BNM and name C11 N1) or (resname DLS1 and name C3 C4 O3 C5 C6 O4 C7 C8 O5 C9 C10 N2 C13 C14 C15) or (resname DLS2 and name C1 C2 O1 C3 C4 O2 C5 C6))) or (segname ${x}T and resname ASP and resid 1 and name CA)"] get index]
     # 4. Define the fixed atom on the flexible chain; no rotations that would move this atom are allowed
     set fa [[atomselect ${molid} "chain ${x} and resname BNM3 and name C18"] get index]
     set tmp [atomselect ${molid} "chain $x"]
#     puts "[$tmp get index]"
#     puts "[$tmp get resname]" 
#     puts "[$tmp get name]"
#     puts "Fixed atom index $fa"
     # Flex-mc in two stages:  1. untargeted; 2. targeted
     # 5. alpha-carbon on the C-terminus of DAVEI Trp3
     set i [[atomselect ${molid} "segname ${x}T and resid [llength $WARHEAD_SEQUENCE] and name CA"] get index]
     puts "i $i"
     # 6. perform an untargeted flex-mc just to remove steric overlaps (see that the $i arg is repeated --
     #    this signals that there is no targeting)
     do_flex_mc $molid $msel $ri $rj $fa $k $i $i $background $rcut $nc $temperature [irand_dom 100 9999] $logid
     # 7. pick alpha-carbon on an env residue near or on MPER
     set pp [atomselect ${molid} "chain $b and resid 648 and name CA"]
     set j [$pp get index]
     # 8. make a phantom target position by moving this atom
     #    attractor_radius Angstroms away from the trimer axis in the XY plane
     set ppr [list [$pp get x] [$pp get y] [$pp get z]]
     set tcom [measure center $background]
     set com [list [lindex $tcom 0] [lindex $tcom 1] [lindex $ppr 2]]
     set dppr [vecsub $ppr $com]
     set rdppr [veclength $dppr]
     set sppr [vecadd $ppr [vecscale $dppr [expr $attractor_radius/$rdppr]]]
     lappend out_targ_points $sppr
     $pp set x [lindex $sppr 0]
     $pp set y [lindex $sppr 1]
     puts "FLEXMC: chain $x end-attracts to index $j at $sppr (r [veclength [vecsub $sppr $com]]) (orig $ppr (r [veclength [vecsub $ppr $com]]))"
     # 9. continue the flex-mc WITH targeting to this phantom atom
     do_flex_mc $molid $msel $ri $rj $fa $k $i $j $background $rcut $nc $temperature [irand_dom 100 9999] $logid
     # 10. restore the phantom CA's x and y positions
     $pp set x [lindex $ppr 0]
     $pp set y [lindex $ppr 1]
     # 11. log the coordinates in the logging molecule
     log_addframe $molid $logid
     $pp delete
    }
    puts "out_targ_points: $out_targ_points"
  }
}


$a writepdb "my_5vn3_mcOut.pdb"

# make a pdb file that fixes all heavy atoms in the original
# crystal structure -- all added atoms are set as unfixed
# for a minimization
mol delete top
mol new my_5vn3.psf
mol addfile my_5vn3_mcOut.pdb
set molid [molinfo top get id]
set a [atomselect ${molid} all]
$a set beta 1
for { set m 0 } { $m < 3 } { incr m } {
   set g [lindex $glist $m]
   set gg [lindex $ggaps $m]
   set ngg [llength $gg]
   for { set s 0 } { $s < $ngg } { incr s } {
      set tg [lindex $gg $s]
      set nrr [lindex $tg 0]
      set crr [lindex $tg 1]
      [atomselect ${molid} "chain $g and resid $nrr to $crr"] set beta 0
   }
   set b [lindex $blist $m]
   set bg [lindex $bgaps $m]
   set nbg [llength $bg]
   for { set s 0 } { $s < $nbg } { incr s } {
      set tg [lindex $bg $s]
      set nrr [lindex $tg 0]
      set crr [lindex $tg 1]
      [atomselect ${molid} "chain $b and resid $nrr to $crr"] set beta 0
   }
}

if { $DOCK_BNM == "1" } {
  [atomselect $molid "same residue as exwithin 4.0 of resname BNM3"] set beta 0
  [atomselect $molid "resname BNM3"] set beta 0
  [atomselect $molid "chain X Y Z"] set beta 0
}

$a writepdb "my_5vn3_fix.pdb"

# finalize logging molecule
# log_finalize $logid
if { $LOG_DCD == "1" } {
   set loga [atomselect $logid all]
   animate write dcd $log_dcd_file waitfor all sel $loga $logid
   mol delete $logid
}

if { [expr $DAVEI > 0] == "1" } {
  # Stage 1 -- SMD-pull trp3's toward the targeting points used above to give clearance
  set tr1 [lindex [[atomselect top "segname XT and resid 1 and name CA"] get {x y z}] 0]
  set tr2 [lindex [[atomselect top "segname YT and resid 1 and name CA"] get {x y z}] 0]
  set tr3 [lindex [[atomselect top "segname ZT and resid 1 and name CA"] get {x y z}] 0]
  set fp [open "clear_daveis_colvars.inp" "w"]
  puts $fp "
  colvarstrajfrequency 100
  colvar {
    name trp3_1
    cartesian {
      atoms {
        psfsegid XT
        atomnameresiduerange CA 9-9
      }
    }
  }
  colvar {
    name trp3_2
    cartesian {
      atoms {
        psfsegid YT
        atomnameresiduerange CA 9-9
      }
    }
  }
  colvar {
    name trp3_3
    cartesian {
      atoms {
        psfsegid ZT
        atomnameresiduerange CA 9-9
      }
    }
  }
  harmonic {
    colvars trp3_1 
    centers ( [lindex $tr1 0] , [lindex $tr1 1] , [lindex $tr1 2] ) 
    targetcenters [format %s [pc3 [lindex $out_targ_points 0]]] 
    forceconstant 1.0
    targetnumsteps 30000 
  }
  harmonic {
    colvars trp3_2 
    centers ( [lindex $tr2 0] , [lindex $tr2 1] , [lindex $tr2 2] ) 
    targetcenters [format %s [pc3 [lindex $out_targ_points 1]]] 
    forceconstant 1.0
    targetnumsteps 30000 
  }
  harmonic {
    colvars trp3_3
    centers ( [lindex $tr3 0] , [lindex $tr3 1] , [lindex $tr3 2] ) 
    targetcenters [format %s [pc3 [lindex $out_targ_points 2]]] 
    forceconstant 1.0
    targetnumsteps 30000 
  }
  "
  close $fp

  # Stage 3: drag Trp3's toward Env MPER's

  if { $MPER_EXTEND == "1" } {
    set env_targ_as "664 to 672"
    set env_targ_cv "664-672"
    set env_avoid_cv "610-630"
  } else {
    set env_targ_as "659 to 663"
    set env_targ_cv "659-663"
    set env_avoid_cv "610-630"
  }
  set dx [veclength [vecsub [lindex $out_targ_points 0] [measure center [atomselect $molid "chain D and resid $env_targ_as"]]]]
  set dy [veclength [vecsub [lindex $out_targ_points 1] [measure center [atomselect $molid "chain A and resid $env_targ_as"]]]]
  set dz [veclength [vecsub [lindex $out_targ_points 2] [measure center [atomselect $molid "chain B and resid $env_targ_as"]]]]
  set fp [open "drag_daveis_colvars.inp" "w"]
  puts $fp "
  colvarstrajfrequency 100
  scriptedcolvarforces on
  colvar {
     name trp3_1
     distance {
       group1 {
         psfsegid XT
         atomnameresiduerange CA 1-9
       }
       group2 {
         psfsegid D
         atomnameresiduerange CA $env_targ_cv
       }
     }
  }
  colvar {
     name trp3_2
     distance {
       group1 {
         psfsegid YT
         atomnameresiduerange CA 1-9
       }
       group2 {
         psfsegid A
         atomnameresiduerange CA $env_targ_cv
       }
     }
  }
  colvar {
     name trp3_3
     distance {
       group1 {
         psfsegid ZT
         atomnameresiduerange CA 1-9
       }
       group2 {
         psfsegid B
         atomnameresiduerange CA $env_targ_cv
       }
     }
  }
  colvar {
     name trp3_1av
     distance {
       group1 {
         psfsegid XT
         atomnameresiduerange CA 1-9
       }
       group2 {
         psfsegid B
         atomnameresiduerange CA $env_avoid_cv
       }
     }
  }
  colvar {
     name trp3_2av
     distance {
       group1 {
         psfsegid YT
         atomnameresiduerange CA 1-9
       }
       group2 {
         psfsegid D
         atomnameresiduerange CA $env_avoid_cv
       }
     }
  }
  colvar {
     name trp3_3av
     distance {
       group1 {
         psfsegid ZT
         atomnameresiduerange CA 1-9
       }
       group2 {
         psfsegid A
         atomnameresiduerange CA $env_avoid_cv
       }
     }
  }
  harmonic {
    colvars trp3_1 trp3_2 trp3_3
    centers $dx $dy $dz
    targetcenters 4.0 4.0 4.0
    forceconstant 10.0
    targetnumsteps 60000 60000 60000
  }
  "
  close $fp
}



# clean up
foreach f $LOCALFILES {
  exec rm $f
}

quit

