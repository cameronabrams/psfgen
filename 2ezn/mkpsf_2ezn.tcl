# VMD/psfgen script for generated psf/pdb pair for PDB 2ezn
# cyanovirin-N
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
for { set i 0 } { $i < [llength $argv] } { incr i } {
  if { [lindex $argv $i] == "-davei" } {
     incr i
     set nSpacerRepeats [lindex $argv $i]
  }
  if { $arg == "-seed" } {
     incr i
     set seed [lindex $argv $i] 
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
package require psfgen
resetpsf
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf

alias residue HIS HSE
alias atom ILE CD1 CD
  
segment A { 
   pdb 2ezn.pdb
   if { [expr $nSpacerRepeats > 0] } { 
      set r 102
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
}
# assign coordinates where available
coordpdb 2ezn.pdb A
# patches!
patch DISU A:8   A:22
patch DISU A:58  A:73

guesscoord

writepsf my_2ezn.psf
writepdb my_2ezn_raw.pdb; lappend LOCALFILES my_2ezn_raw.pdb

mol new my_2ezn.psf 
mol addfile my_2ezn_raw.pdb
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
$a writepdb my_2ezn.pdb

$a set beta 0
set b [atomselect top "name CA"]
$b set beta 1
$a writepdb my_2ezn_fix.pdb

quit
