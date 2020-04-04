# VMD/psfgen script for generating psf/pdb pair for PDB 6w41
# SARS-CoV-2 S RBD and CR3022 Fab
#
# cameron f abrams (c) 2020
# drexel university
# chemical and biological engineering
#
# check for base directory name variable;
# if not set, use default

set BASEPDB 6w41

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
  if  { $arg == "-rbdonly" } {
    set RBDONLY 1
  }
}

expr srand($seed)

# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

set DOMC 0

mol new ${BASEPDB}.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

##### output of python3 parse_pdb_psfgen.py -pdb 6w41.pdb below #####a
[atomselect top "chain C and protein and resid 333 to 527"] writepdb "C_333_to_527.pdb"
segment C {
   pdb C_333_to_527.pdb
}
coordpdb C_333_to_527.pdb C
if { $RBDONLY == 0 } {
[atomselect top "chain H and protein and resid 1 to 216"] writepdb "H_1_to_216.pdb"
segment H {
   pdb H_1_to_216.pdb
}
coordpdb H_1_to_216.pdb H
[atomselect top "chain L and protein and resid 1 to 215"] writepdb "L_1_to_215.pdb"
segment L {
   pdb L_1_to_215.pdb
}
coordpdb L_1_to_215.pdb L
}
set myseg [atomselect top "chain C and resid 601 to 601"]
$myseg set resname BGNA
$myseg writepdb CS_601_to_601.pdb
segment CS {
   pdb CS_601_to_601.pdb
}
coordpdb CS_601_to_601.pdb CS
if { $RBDONLY == 0} {
patch DISU H:22 H:92
patch DISU H:140 H:196
patch DISU H:216 L:214
patch DISU L:23 L:88
patch DISU L:134 L:194
}
patch DISU C:336 C:361
patch DISU C:379 C:432
patch DISU C:391 C:525
patch DISU C:480 C:488
patch NGLB C:343 CS:601
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

