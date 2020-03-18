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

#pdbalias residue NAG BGLC
#pdbalias atom BGLC C7 C
#pdbalias atom BGLC O7 O
#pdbalias atom BGLC C8 CT

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
set I [atomselect top "ions"] 
$I set resname ZN2
$I writepdb "ions.pdb"
segment I {
   pdb ions.pdb
}
coordpdb ions.pdb I
patch DISU A:133 A:141
patch DISU A:344 A:361
patch DISU A:530 A:542
patch DISU E:336 E:361
patch DISU E:379 E:432
patch DISU E:391 E:525
patch DISU E:480 E:488
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

