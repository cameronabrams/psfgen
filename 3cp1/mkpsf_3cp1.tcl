# VMD/psfgen script for generating psf/pdb pair for PDB 3cp1
# HIV-1 gp41 NHR/CHR six-helix bundle with DL replaced with GGG
#
# cameron f abrams (c) 2017-2019
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
set LOG_DCD 0
set logid -1
# NHR-CHR linker sequence
set NHR_CHR_linker [list GLY GLY GLY GLY GLY]
set seed 12345
# if set, this flag indicates that residues 663 and beyond are to be modeled in
# (unimplemented)
set MPER_EXTEND 0
set TMD_EXTEND 0
set MPER_663_to_683 [list LEU ASP LYS TRP ALA SER LEU TRP ASN TRP PHE ASN ILE THR ASN TRP LEU TRP TYR ILE LYS]
set TMD_684_to_709 [list LEU PHE ILE MET ILE VAL GLY GLY LEU VAL GLY LEU ARG ILE VAL PHE ALA VAL LEU SER ILE VAL ASN ARG VAL ARG]
set SKIP_LOOPMC 0
set DO_STALK 0
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
  if { $arg == "-tmd-extend" } {
     set MPER_EXTEND 1
     set TMD_EXTEND 1
  }
   if { $arg == "-skip-loopmc" } {
     set SKIP_LOOPMC 1
  }
  if { $arg == "-log-dcd" } {
     set LOG_DCD 1
     incr a
     set log_dcd_file [lindex $argv $a]
  }
  if { $arg == "-do-stalk" } {
     set DO_STALK 1
  }
}

expr srand($seed)

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

mol new 3cp1.pdb
set a [atomselect top all]
# renumber to HXB2
set n [atomselect top "protein and resid 3 to 42"]
set rid [$n get resid]
for {set i 0} {$i < [llength $rid]} {incr i} {
   lset rid $i [expr [lindex $rid $i] + 535]
}
$n set resid $rid
$n set chain A
$n set segname A
$n writepdb "NHR-A.pdb"; lappend LOCALFILES NHR-A.pdb
set c [atomselect top "protein and resid 52 to 85"]
set rid [$c get resid]
for {set i 0} {$i < [llength $rid]} {incr i} {
   lset rid $i [expr [lindex $rid $i] + 577]
}
$c set resid $rid
$c set chain A
$c set segname A
$c writepdb "CHR-A.pdb"; lappend LOCALFILES CHR-A.pdb

# get the waters that fully occupy cell
set w1 [atomselect top "resname HOH and occupancy 1.00"]
$w1 set name OH2
$w1 set resname TIP3
$w1 set segname WX1
$w1 set chain WX1
$w1 writepdb "WX1.pdb"; lappend LOCALFILES WX1.pdb

# get the waters that 1/3 occupy cell (these will not be replicated by symmetry in
# building the trimer)
set w0 [atomselect top "resname HOH and occupancy 0.33"]
$w0 set name OH2
$w0 set resname TIP3
$w0 set segname WX0
$w0 set chain WX0
$w0 writepdb "WX0.pdb"; lappend LOCALFILES WX0.pdb

# save position of backbone N that will define placement of first modeled-in residue
set n577pos(A) [cacoIn_nOut 576 A 0]
set n663pos(A) [cacoIn_nOut 662 A 0]

# BIOMT 2 from pdb
set tmat3 {{0 1 0 0} {0 0 1 0} {1 0 0 0} {0 0 0 1}}
foreach ch { B C } N { 2 3 } {
   $a move $tmat3
   $n set chain $ch
   $n set segname $ch
   $c set chain $ch
   $c set segname $ch
   $w1 set chain WX${N}
   $w1 set segname WX${N}
   $n writepdb "NHR-${ch}.pdb";lappend LOCALFILES NHR-${ch}.pdb
   $c writepdb "CHR-${ch}.pdb";lappend LOCALFILES CHR-${ch}.pdb
   set n577pos($ch) [cacoIn_nOut 576 $ch 0]
   set n663pos($ch) [cacoIn_nOut 662 $ch 0]
   $w1 writepdb "WX${N}.pdb";lappend LOCALFILES WX${N}.pdb
}

mol delete top
package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str
pdbalias residue HIS HSD
resetpsf

foreach s {A B C} {
  segment ${s} {
    pdb NHR-${s}.pdb
    set r 577
    foreach aa $NHR_CHR_linker {
	residue $r $aa $s
	incr r
    }
    set NHR_CHR_linker_Cterm [expr $r - 1]
    pdb CHR-${s}.pdb
    if { $MPER_EXTEND == "1" } {
       set lr 0
       for { set r 663 } { $r < 684 } { incr r } {
          residue $r [lindex $MPER_663_to_683 $lr] ${s}
	  incr lr
       }
    } 
    if { $TMD_EXTEND == "1" } { 
       set lr 0
       for { set r 684 } { $r < 710 } { incr r } {
	  puts "[lindex $TMD_684_to_709 $lr]"
          residue $r [lindex $TMD_684_to_709 $lr] $s
	  incr lr
       }
    }
  }
  coordpdb NHR-${s}.pdb ${s}
  coordpdb CHR-${s}.pdb ${s}
  coord ${s} 577 N $n577pos(${s})
  if { $MPER_EXTEND } {
    coord ${s} 663 N $n663pos(${s})
  }
}

foreach w {0 1 2 3} {
  segment WX${w} {
     auto none
     pdb WX${w}.pdb
  }
  coordpdb WX${w}.pdb WX${w}
}

guesscoord
regenerate angles dihedrals
writepsf "my_3cp1.psf"
writepdb "my_3cp1_stage1.pdb"
lappend LOCALFILES my_3cp1_stage1.pdb

mol new my_3cp1.psf
mol addfile my_3cp1_stage1.pdb
set molid [molinfo top get id]

# rotate to get trimer axis along z
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
set ca [measure center [atomselect top "protein and resid 538 to 546"] weight mass]
set cb [measure center [atomselect top "protein and resid 568 to 676"] weight mass]
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
$a writepdb my_3cp1_stage2.pdb

# initialize a logging molecule
if { $LOG_DCD == "1" } {
   # set logid [log_initialize my_${SYSNAME}.psf my_${SYSNAME}.pdb]
   mol new my_3cp1.psf
   mol addfile my_3cp1_stage2.pdb
   set logid [molinfo top get id]
   puts "Logging molecule number $logid to be saved to $log_dcd_file" 
   mol top $molid
}
set a [atomselect top all]
$a set beta 1 ;# pretty much everything is fixed except...

# fold the optional MPER and TM extensions into alpha helices
if { $MPER_EXTEND == "1" } {
  foreach c {A B C} {
    set Cterm 683
    if { $TMD_EXTEND == "1" } {
      set Cterm 709
    }
    set sel [atomselect $molid "protein and chain $c and resid 661 to $Cterm"]
    puts "Fold-alpha: chain $c MPER/TM extension"
    fold_alpha_helix $molid $sel extrabonds-${c}.inp
    lappend LOCALFILES extrabonds-${c}.inp
    $sel delete
    log_addframe $molid $logid
    
    if { $DO_STALK == "1" } {
       # fold them over!
       puts "Folding 'em over!..."
       set endres [[atomselect $molid "protein and chain $c and resid $Cterm and name CA"] get residue]
       set rot1res [[atomselect $molid "protein and chain $c and resid 662 and name CA"] get residue]
       set rot2res [[atomselect $molid "protein and chain $c and resid 672 and name CA"] get residue]
       Crot_psi_toCterm $rot1res $endres $c $molid 120
       log_addframe $molid $logid
       Crot_phi_toCterm $rot1res $endres $c $molid -60
       log_addframe $molid $logid
       Crot_psi_toCterm $rot2res $endres $c $molid 120
       log_addframe $molid $logid
    }
  }
}

$a writepdb "my_3cp1_mcIn.pdb"
lappend LOCALFILES my_3cp1_mcIn.pdb

if { $SKIP_LOOPMC == "0" } {
   set nc $LOOP_MC_NC
   set rcut $LOOP_MC_RCUT
   set temperature $LOOP_MC_TEMPERATURE
   set k $LOOP_MC_K
   set r0 $LOOP_MC_R0
   set background [atomselect ${molid} "noh"]
   foreach ch { A B C } {
     puts "MC on chain $ch of molid $molid..."
     set residueList [[atomselect ${molid} "chain ${ch} and resid 577 to $NHR_CHR_linker_Cterm and name CA"] get residue]
     do_loop_mc ${residueList} ${ch} ${molid} ${k} ${r0} ${background} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] $logid 
   }
}

$a writepdb "my_3cp1_mcOut.pdb"
if { $DO_STALK == "1" } {
  set fixem [atomselect top "(protein and noh and not resid 577 to 579 661 662 663 671 672 673) or (not protein and name OH2)"]
} else {
  set fixem [atomselect top "(protein and noh and not resid 577 to 579) or (not protein and name OH2)"]
}
  $a set beta 0
$fixem set beta 1
$a writepdb "my_3cp1_fixed.pdb"


# finalize logging molecule
# log_finalize $logid
if { $LOG_DCD == "1" } {
   set loga [atomselect $logid all]
   animate write dcd $log_dcd_file waitfor all sel $loga $logid
   mol delete $logid
}

foreach f $LOCALFILES {
  exec rm $f
}

quit

