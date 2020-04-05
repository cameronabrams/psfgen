# VMD/psfgen script for generating psf/pdb pair for PDB 6waq
# SARS-CoV-2 S RBD and nanobody VHH-72
#
# cameron f abrams (c) 2020
# drexel university
# chemical and biological engineering
#
# check for base directory name variable;
# if not set, use default

set BASEPDB 6waq

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
set RDBONLY 0
set chains { A B C D }
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
  if  { $arg == "-rbdonlyA" } {
    set RBDONLY 1
    set chains { A }
  }
  if  { $arg == "-rbdonlyC" } {
    set RBDONLY 1
    set chains { C }
  }
  if { $arg == "-chainsAB" } {
    set chains { A B }
  }
  if { $arg == "-chainsCD" } {
    set chains { C D }
  }
}

expr srand($seed)

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

mol new ${BASEPDB}.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

##### output of python3 parse_pdb_psfgen.py -pdb 6w41.pdb below #####a
if { [lsearch $chains "A"] != -1 } {
[atomselect top "chain A and protein and resid 1 to 115"] writepdb "A_1_to_115.pdb"
segment A {
   pdb A_1_to_115.pdb
}
coordpdb A_1_to_115.pdb A
}
if { [lsearch $chains "B"] != -1 } {
[atomselect top "chain B and protein and resid 321 to 510"] writepdb "B_321_to_510.pdb"
segment B {
   pdb B_321_to_510.pdb
}
coordpdb B_321_to_510.pdb B
}
if { [lsearch $chains "C"] != -1 } {
[atomselect top "chain C and protein and resid 1 to 113"] writepdb "C_1_to_113.pdb"
segment C {
   pdb C_1_to_113.pdb
}
coordpdb C_1_to_113.pdb C
}
if { [lsearch $chains "D"] != -1 } {
[atomselect top "chain D and protein and resid 320 to 503"] writepdb "D_320_to_503.pdb"
segment D {
   pdb D_320_to_503.pdb
}
coordpdb D_320_to_503.pdb D
set myseg [atomselect top "chain D and resid 601 to 602"]
$myseg set resname BGNA
$myseg writepdb DS_601_to_602.pdb
segment DS {
   pdb DS_601_to_602.pdb
}
coordpdb DS_601_to_602.pdb DS
}
if { [lsearch $chains "B"] != -1 } {
set myseg [atomselect top "chain B and resid 601 to 601"]
$myseg set resname BGNA
$myseg writepdb BS_601_to_601.pdb
segment BS {
   pdb BS_601_to_601.pdb
}
coordpdb BS_601_to_601.pdb BS
}
if { [lsearch $chains "A"] != -1 } {
set myseg [atomselect top "chain A and resid 201 to 203"]
$myseg set name OH2
$myseg set resname TIP3
$myseg writepdb A_201_to_203.pdb
segment AWX {
   pdb A_201_to_203.pdb
}
coordpdb A_201_to_203.pdb AWX
}
if { [lsearch $chains "D"] != -1 }  {
set myseg [atomselect top "chain D and resid 701 to 750"]
$myseg set name OH2
$myseg set resname TIP3
$myseg writepdb D_701_to_750.pdb
segment DWX {
   pdb D_701_to_750.pdb
}
coordpdb D_701_to_750.pdb DWX
set myseg [atomselect top "chain C and resid 201 to 227"]
$myseg set name OH2
$myseg set resname TIP3
$myseg writepdb C_201_to_227.pdb
segment CWX {
   pdb C_201_to_227.pdb
}
coordpdb C_201_to_227.pdb CWX
}
if { [lsearch $chains "B"] != -1 } {
set myseg [atomselect top "chain B and resid 701 to 722"]
$myseg set name OH2
$myseg set resname TIP3
$myseg writepdb B_701_to_722.pdb
segment BWX {
   pdb B_701_to_722.pdb
}
coordpdb B_701_to_722.pdb BWX
}

if { [lsearch $chains "A"] != -1 } {
patch DISU A:22 A:92
}
if { [lsearch $chains "D"] != -1 } {
patch DISU D:323 D:348
patch DISU D:366 D:419
patch DISU D:467 D:474
}
if { [lsearch $chains "C"] != -1 } {
patch DISU C:22 C:92
}
if { [lsearch $chains "B"] != -1 } {
patch DISU B:323 B:348
patch DISU B:366 B:419
patch DISU B:467 B:474
}
if { [lsearch $chains "B"] != -1 && [lsearch $chains "D"] != -1 } {
patch DISU D:378 B:378
}
if { [lsearch $chains "D"] != -1 } {
patch NGLB D:330 DS:601
patch NGLB D:357 DS:602
}
if { [lsearch $chains "B"] != -1 } {
patch NGLB B:330 BS:601
}
##### output of python3 parse_pdb_psfgen.py -pdb 6w41.pdb above #####

guesscoord

regenerate angles dihedrals

writepsf "my_${BASEPDB}.psf"
writepdb "unrelaxed.pdb"


mol delete top
mol new my_${BASEPDB}.psf
mol addfile unrelaxed.pdb
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
$a writepdb "my_${BASEPDB}.pdb"

# clean up
foreach f $LOCALFILES { 
  exec rm $f
}

quit

