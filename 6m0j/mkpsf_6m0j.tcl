# VMD/psfgen script for generating psf/pdb pair for PDB 6m0j
# complex of SARS-CoV-2 S RBD and ACE2
#
# cameron f abrams (c) 2020
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

set seed 12345
set LOG_DCD 0
set logid -1
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-seed" } {
    incr a
    set seed [lindex $argv $a]
  }
  if { $arg == "-log-dcd" } {
    set LOG_DCD 1
    incr a
    set log_dcd_file [lindex $argv $a]
  }
}

expr srand($seed)

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

set DOMC 0

mol new 6m0j.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

pdbalias residue NAG BGNA
pdbalias atom BGNA C7 C
pdbalias atom BGNA O7 O
pdbalias atom BGNA C8 CT

##### output of python3 parse_pdb_psfgen.py 6m0j.pdb below #####a
[atomselect top "chain A and protein and resid 19 to 615"] writepdb "A_19_to_615.pdb"
segment A {
   pdb A_19_to_615.pdb
}
coordpdb A_19_to_615.pdb A
[atomselect top "chain E and protein and resid 333 to 526"] writepdb "E_333_to_526.pdb"
segment E {
   pdb E_333_to_526.pdb
}
coordpdb E_333_to_526.pdb E
set myseg [atomselect top "chain A and resid 901 to 901"]
$myseg set resname ZN2
$myseg writepdb "AI_901_to_901.pdb"
segment AI {
   pdb AI_901_to_901.pdb
}
coordpdb AI_901_to_901.pdb AI
set myseg [atomselect top "chain A and resid 902 to 904"]
$myseg set resname BGNA
$myseg writepdb AS_902_to_904.pdb
segment AS {
   pdb AS_902_to_904.pdb
}
coordpdb AS_902_to_904.pdb AS
set myseg [atomselect top "chain E and resid 601 to 601"]
$myseg set resname BGNA
$myseg writepdb ES_601_to_601.pdb
segment ES {
   pdb ES_601_to_601.pdb
}
coordpdb ES_601_to_601.pdb ES
set myseg [atomselect top "chain A and resid 1001 to 1066"]
$myseg set name OH2
$myseg set resname TIP3
$myseg writepdb A_1001_to_1066.pdb
segment AWX {
   pdb A_1001_to_1066.pdb
}
coordpdb A_1001_to_1066.pdb AWX
set myseg [atomselect top "chain E and resid 701 to 714"]
$myseg set name OH2
$myseg set resname TIP3
$myseg writepdb E_701_to_714.pdb
segment EWX {
   pdb E_701_to_714.pdb
}
coordpdb E_701_to_714.pdb EWX
patch DISU A:133 A:141
patch DISU A:344 A:361
patch DISU A:530 A:542
patch DISU E:336 E:361
patch DISU E:379 E:432
patch DISU E:391 E:525
patch DISU E:480 E:488
patch NGLB A:90 AS:902
patch NGLB A:322 AS:904
patch NGLB A:546 AS:903
patch NGLB E:343 ES:601
##### output of python3 parse_pdb_psfgen.py 6m0j.pdb above #####

guesscoord

regenerate angles dihedrals

writepsf "my_6m0j.psf"
writepdb "unrelaxed.pdb"

mol delete top
mol new my_6m0j.psf
mol addfile unrelaxed.pdb
set molid [molinfo top get id]
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
$a writepdb "my_6m0j.pdb"

# clean up
foreach f $LOCALFILES { 
  exec rm $f
}

quit

