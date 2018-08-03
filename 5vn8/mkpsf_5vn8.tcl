# VMD/psfgen script for generating psf/pdb pair for PDB 5vn8
# trimeric HIV-1 gp140 in "closed" state
# 
# note: b12 FAbs antibodies (6 chains) are not included as yet.
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

set MPER_665_to_682 [list LYS TRP ALA SER LEU TRP ASN TRP PHE ASP ILE SER ASN TRP LEU TRP TYR ILE]
set MPER_EXTEND 0
set seed 12345
set SKIP_LOOPMC 0
set LOG_DCD 0
set logid -1
set DOCK_BNM 0
set LOOP_MC_NC  1000
set LOOP_MC_RCUT  3.0
set LOOP_MC_TEMPERATURE  2.5
set LOOP_MC_K  10.0
set LOOP_MC_R0 1.5 
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
}

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

set SYSNAME 5vn8
expr srand($seed)

mol new 5vn8.pdb
set molid [molinfo top get id]

# extract contiguous protein segments as individual pdb's                                          
set glist { G E D }
set blist { A B C }
set xlist { X Y Z }
set gsegs {
  {{32 58} {70 117} {123 136} {149 165} {175 204} {209 505}} 
  {{32 58} {70 117} {123 136} {149 165} {175 204} {209 505}} 
  {{32 58} {70 117} {123 136} {149 165} {175 204} {209 505}} }
set bsegs {
  {{512 549} {564 664}} 
  {{512 549} {564 664}} 
  {{512 549} {564 664}} }

# identify missing residues
set ggaps {
  {{  59   69 { LYS ALA TYR ASP THR GLU VAL HSD ASN VAL TRP }} 
   { 118  122 { PRO CYS VAL LYS LEU }} 
   {136A 136P { ASN ASN THR ASN ASN SER THR ASN ALA THR ILE SER ASP TRP GLU LYS }}
   { 166  174 { ARG ASP LYS ILE LYS LYS GLU TYR ALA }}}
  {{  59   69 { LYS ALA TYR ASP THR GLU VAL HSD ASN VAL TRP }} 
   { 118  122 { PRO CYS VAL LYS LEU }} 
   {136A 136P { ASN ASN THR ASN ASN SER THR ASN ALA THR ILE SER ASP TRP GLU LYS }}
   { 166  174 { ARG ASP LYS ILE LYS LYS GLU TYR ALA }}}
  {{  59   69 { LYS ALA TYR ASP THR GLU VAL HSD ASN VAL TRP }} 
   { 118  122 { PRO CYS VAL LYS LEU }} 
   {136A 136P { ASN ASN THR ASN ASN SER THR ASN ALA THR ILE SER ASP TRP GLU LYS }}
   { 166  174 { ARG ASP LYS ILE LYS LYS GLU TYR ALA }}}}

set bgaps {
  {{ 550  563 { GLN GLN GLN ASN ASN LEU LEU ARG ALA PRO GLU ALA GLN GLN }}}
  {{ 550  563 { GLN GLN GLN ASN ASN LEU LEU ARG ALA PRO GLU ALA GLN GLN }}}
  {{ 550  563 { GLN GLN GLN ASN ASN LEU LEU ARG ALA PRO GLU ALA GLN GLN }}}}

# identify disulfide
set ssbonds {
  { G  54 G  74 } { G 126 G 196 } { G 131 G 157 } { G 218 G 247 } { G 228 G 239 } 
  { G 296 G 331 } { G 378 G 445 } { G 385 G 418 } { G 501 A 605 } { A 598 A 604 } 
  { E  54 E  74 } { E 126 E 196 } { E 131 E 157 } { E 218 E 247 } { E 228 E 239 } 
  { E 296 E 331 } { E 378 E 445 } { E 385 E 418 } { E 501 C 605 } { C 598 C 604 } 
  { D  54 D  74 } { D 126 D 196 } { D 131 D 157 } { D 218 D 247 } { D 228 D 239 } 
  { D 296 D 331 } { D 378 D 445 } { D 385 D 418 } { D 501 B 605 } { B 598 B 604 } 
}

set mln {}
foreach g $glist b $blist {
  lappend mln [list $g 59 58 ]
  lappend mln [list $g 118 117 ]
  lappend mln [list $g 136A 136 ]
  lappend mln [list $g 166 165 ]
  lappend mln [list $b 550 549 ] 
}

set ns {}
foreach m $mln {
  lappend ns [cacoIn_nOut [lindex $m 2] [lindex $m 0] 0]
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
foreach chain { G E D A B C } {
  set g [atomselect $molid "chain $chain and not protein"]
  [atomselect $molid "chain $chain and resname NAG"] set resname BGNA
  [atomselect $molid "chain $chain and resname BMA"] set resname BMAN
  [atomselect $molid "chain $chain and resname MAN"] set resname AMAN
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

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias atom BGNA N2 N
pdbalias atom BGNA C7 C
pdbalias atom BGNA O7 O
pdbalias atom BGNA C8 CT
foreach nbad [exec grep 5VG 5f4p.pdb | grep -w "5VG A" | grep HETATM | cut -b 13-16] ngood [exec grep ATOM ${PSFGEN_BASEDIR}/charmm/bnm.str | cut -b 6-9 | grep -v ^H] {
  pdbalias atom BNM3 $nbad $ngood
}

proc my_increment { numlet } {
   if { [string is integer $numlet] == "1" } {
     return [expr $numlet + 1]
   }
   set c [format %c [expr [scan [string index $numlet end] %c] + 1]]
   set nn [string replace $numlet end end $c]
   return $nn
}

for { set m 0 } { $m < 3 } { incr m } {
  set g [lindex $glist $m]
  set b [lindex $blist $m]

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
    set x [lindex $xlist $m]
    segment ${x} {
      pdb ${x}.pdb
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

# set N positions for first residue of each modeled-in loop
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
  patch NGLB ${g}:88 ${g}S:601
  patch NGLB ${g}:160 ${g}S:602
  patch NGLB ${g}:197 ${g}S:603
  patch NGLB ${g}:234 ${g}S:605
  patch NGLB ${g}:241 ${g}S:607
  patch NGLB ${g}:262 ${g}S:610
  patch NGLB ${g}:276 ${g}S:614
  patch NGLB ${g}:295 ${g}S:616
  patch NGLB ${g}:301 ${g}S:618
  patch NGLB ${g}:332 ${g}S:619
  patch NGLB ${g}:339 ${g}S:621
  patch NGLB ${g}:355 ${g}S:622
  patch NGLB ${g}:362 ${g}S:623
  patch NGLB ${g}:386 ${g}S:626
  patch NGLB ${g}:392 ${g}S:629
  patch NGLB ${g}:397 ${g}S:630
  patch NGLB ${g}:413 ${g}S:631
  patch NGLB ${g}:448 ${g}S:633

  patch NGLB ${b}:611 ${b}S:701
  patch NGLB ${b}:616 ${b}S:702
  patch NGLB ${b}:625 ${b}S:703
  patch NGLB ${b}:637 ${b}S:704

  patch 13ab ${g}S:612 ${g}S:613

  patch 14bb ${g}S:603 ${g}S:604
  patch 14bb ${g}S:605 ${g}S:606
  patch 14bb ${g}S:607 ${g}S:608
  patch 14bb ${g}S:608 ${g}S:609
  patch 14bb ${g}S:610 ${g}S:611
  patch 14bb ${g}S:611 ${g}S:612
  patch 14bb ${g}S:614 ${g}S:615
  patch 14bb ${g}S:616 ${g}S:617
  patch 14bb ${g}S:619 ${g}S:620
  patch 14bb ${g}S:623 ${g}S:624
  patch 14bb ${g}S:624 ${g}S:625
  patch 14bb ${g}S:626 ${g}S:627
  patch 14bb ${g}S:627 ${g}S:628
  patch 14bb ${g}S:631 ${g}S:632
  patch 14bb ${g}S:633 ${g}S:634

}

guesscoord
regenerate angles dihedrals

writepsf "my_5vn8.psf"
writepdb "unrelaxed.pdb"

resetpsf

lappend LOCALFILES unrelaxed.pdb

mol new my_5vn8.psf
mol addfile unrelaxed.pdb

set molid [molinfo top get id]

# initialize a logging molecule
if { $LOG_DCD == "1" } {
   # set logid [log_initialize my_${SYSNAME}.psf my_${SYSNAME}.pdb]
   mol new my_${SYSNAME}.psf
   mol addfile unrelaxed.pdb
   set logid [molinfo top get id]
   puts "Logging molecule number $logid to be saved to $log_dcd_file"
   mol top $molid
}

if { $MPER_EXTEND == "1" } {
   for { set m 0 } { $m < 3 } { incr m } {
     set b [lindex $blist $m]
     set Cterm 682
     set sel [atomselect $molid "protein and chain $b and resid 664 to $Cterm"]
     puts "Fold-alpha: chain $b MPER extension"
     fold_alpha_helix $molid $sel extrabonds-${b}.inp
     log_addframe $molid $logid
     $sel delete
   }
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
       do_loop_mc ${residueList} $b ${molid} ${k} ${r0} ${background} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] $logid
    }
  }
}

set a [atomselect ${molid} all]

$a writepdb "my_5vn8_mcOut.pdb"

# make a pdb file that fixes all heavy atoms in the original
# crystal structure -- all added atoms are set as unfixed
# for a minimization
mol delete top
mol new my_5vn8.psf
mol addfile my_5vn8_mcOut.pdb
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
}

$a writepdb "my_5vn8_fix.pdb"

# finalize logging molecule
# log_finalize $logid
if { $LOG_DCD == "1" } {
   set loga [atomselect $logid all]
   animate write dcd $log_dcd_file waitfor all sel $loga $logid
   mol delete $logid
}

# clean up
foreach f $LOCALFILES {
  exec rm $f
}

quit

