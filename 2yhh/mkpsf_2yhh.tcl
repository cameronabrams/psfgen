# VMD/psfgen script for generated psf/pdb pair for PDB 2yhh
# microvirin
#
# option to append a spacer, his-tag, and MPER sequence to the C-terminus
# to make a CVN-Linker-MPER DAVEI.
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

# check argument list
set nSpacerRepeats 0
set seed 12445
set mutate 0
set trp3
for { set i 0 } { $i < [llength $argv] } { incr i } {
  if { [lindex $argv $i] == "-davei" } {
     incr i
     set nSpacerRepeats [lindex $argv $i]
  }
  if { [lindex $argv $i] == "-trp3" } {
     set trp3 1
  }
  if { [lindex $argv $i] == "-seed" } {
     incr i
     set seed [lindex $argv $i] 
  }
  if { [lindex $argv $i] == "-mutate" } {
     set mutate 1
  }
}

expr srand($seed)

# Make fusion with MPER c-terminal to CVN
set SPACER_SEQ { }
for { set i 0 } { $i < $nSpacerRepeats } { incr i } {
  lappend SPACER_SEQ "GLY GLY GLY SER"
}
set SPACER_SEQ [join $SPACER_SEQ]
puts "SPACER: $SPACER_SEQ"
set H6 {HSE HSE HSE HSE HSE HSE}
set MPER_SEQ { ASP LYS TRP ALA SER LEU TRP ASN TRP PHE GLU ILE THR GLU TRP LEU TRP TYR ILE LYS }
if { $trp3 == "1" } {
  set MPER_SEQ { ASP LYS TRP ALA SER LEU TRP ASN TRP }
}

mol new 2yhh.pdb
set a [atomselect top "protein"]
$a writepdb "2yhh_prot.pdb"
set a [atomselect top "resname MAN"]
$a set resname AMAN
$a writepdb "2yhh_man.pdb"

package require psfgen
resetpsf
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology /home/cfa/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology /home/cfa/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str

alias residue HIS HSE
alias atom ILE CD1 CD

segment A { 
   pdb 2yhh_prot.pdb
   if { [expr $nSpacerRepeats > 0] } { 
      set r 109
      set R0 $r
      for { set i 0 } { $i < [llength $SPACER_SEQ] } { incr i } {
        residue [expr $r + $i] [lindex $SPACER_SEQ $i] A
      }
      set r [expr $r + [llength $SPACER_SEQ]]
      for { set i 0 } { $i < [llength $H6] } { incr i } {
        residue [expr $r + $i] [lindex $H6 $i] A
      }
      set r [expr $r + [llength $H6]]	
      set R1 $r
      for { set i 0 } { $i < [llength $MPER_SEQ] } { incr i } {
        residue [expr $r + $i] [lindex $MPER_SEQ $i] A
      }
      set R2 [expr $r + $i]
   }
   if { $mutate } {
     mutate 81 LYS
     mutate 83 ARG
   }
}
segment M {
   pdb 2yhh_man.pdb
}

# assign coordinates where available
coordpdb 2yhh_prot.pdb A
coordpdb 2yhh_man.pdb M
# patches!
patch DISU A:8   A:24
patch DISU A:60  A:80
patch DISU A:63  A:78
patch 12aa M:1109 M:1110

regenerate angles dihedrals

guesscoord

writepsf my_2yhh.psf
writepdb my_2yhh_raw.pdb; lappend LOCALFILES my_2yhh_raw.pdb

mol new my_2yhh.psf 
mol addfile my_2yhh_raw.pdb
set a [atomselect top all]
if { [expr $nSpacerRepeats > 0] } {
  set ri [[atomselect top "resid $R0 to $R2 and name N CA"] get index]
  set rj [[atomselect top "resid $R0 to $R2 and name CA C"] get index]
  set fa [[atomselect top "resid 100 and name CA"] get index]
  set bg [atomselect top "resid $R0 to $R2"]
  do_flex_mc top $a $ri $rj $fa 0 0 0 $bg 3.0 1000 2.5 121
  set b [atomselect top "resid $R1 to $R2"]
  fold_alpha_helix top $b
}
$a writepdb my_2yhh.pdb

$a set beta 0
set b [atomselect top "name CA"]
$b set beta 1
$a writepdb my_2yhh_fix.pdb

quit
