# VMD/psfgen script for generated psf/pdb pair for PDB 5jyn
# HIV-1 gp41 TM trimer in DMPC bilayer
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

# check argument list
set seed 12445
for { set i 0 } { $i < [llength $argv] } { incr i } {
  if { [lindex $argv $i] == "-seed" } {
     incr i
     set seed [lindex $argv $i] 
  }
}

expr srand($seed)
set LOCALFILES {}

mol new 5jyn.pdb waitfor all
set a [atomselect top all]
$a frame 14
$a writepdb "5jyn_1.pdb"; lappend LOCALFILES 5jyn_1.pdb
mol delete top
mol new 5jyn_1.pdb
alnpa top
set a [atomselect top all] 
$a move [transaxis x 180 degrees]
$a moveby [vecscale -1 [measure center $a]]

foreach c {A B C} {
  [atomselect top "chain $c"] writepdb "${c}.pdb"; lappend LOCALFILES ${c}.pdb
}
mol delete top

package require psfgen
resetpsf
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
alias residue HIS HSE
alias atom ILE CD1 CD
  
foreach s {A B C} {
  segment ${s} { 
    pdb ${s}.pdb
  }
}
foreach s {A B C} {
  coordpdb ${s}.pdb $s
}

guesscoord

writepsf my_5jyn.psf
writepdb my_5jyn.pdb

foreach f $LOCALFILES {
  exec rm -f $f
}

quit
