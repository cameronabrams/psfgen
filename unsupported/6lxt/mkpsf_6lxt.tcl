# VMD/psfgen script for generating psf/pdb pair for PDB 6lxt
# SARS-CoV-2 S2 6hb construct
#
# cameron f abrams (c) 2020
# drexel university
# chemical and biological engineering
#
# check for base directory name variable;
# if not set, use default

set BASEPDB 6lxt

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

mol new ${BASEPDB}.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

##### output of python3 parse_pdb_psfgen.py 6m0j.pdb below #####a
set segs  { { A  912  988 }  { A 1164 1202 } 
            { B  914  988 }  { B 1162 1202 } 
            { C  912  988 }  { C 1163 1202 }  } 
set loops { { A 1156 1163 } 
            { B 1156 1161 } 
            { C 1156 1162 }  }
[atomselect top "chain A and protein and resid 912 to 988"] writepdb "A_912_to_988.pdb"
[atomselect top "chain A and protein and resid 1164 to 1202"] writepdb "A_1164_to_1202.pdb"
segment A {
   pdb A_912_to_988.pdb
   residue 1156 SER A
   residue 1157 GLY A
   residue 1158 GLY A
   residue 1159 ARG A
   residue 1160 GLY A
   residue 1161 GLY A
   residue 1162 PRO A
   residue 1163 ASP A
   pdb A_1164_to_1202.pdb
}
coordpdb A_912_to_988.pdb A
coordpdb A_1164_to_1202.pdb A
coord A 1156 N [cacoIn_nOut 988 A 0]
[atomselect top "chain B and protein and resid 914 to 988"] writepdb "B_914_to_988.pdb"
[atomselect top "chain B and protein and resid 1162 to 1202"] writepdb "B_1162_to_1202.pdb"
segment B {
   pdb B_914_to_988.pdb
   residue 1156 SER B
   residue 1157 GLY B
   residue 1158 GLY B
   residue 1159 ARG B
   residue 1160 GLY B
   residue 1161 GLY B
   pdb B_1162_to_1202.pdb
}
coordpdb B_914_to_988.pdb B
coordpdb B_1162_to_1202.pdb B
coord B 1156 N [cacoIn_nOut 988 B 0]
[atomselect top "chain C and protein and resid 912 to 988"] writepdb "C_912_to_988.pdb"
[atomselect top "chain C and protein and resid 1163 to 1202"] writepdb "C_1163_to_1202.pdb"
segment C {
   pdb C_912_to_988.pdb
   residue 1156 SER C
   residue 1157 GLY C
   residue 1158 GLY C
   residue 1159 ARG C
   residue 1160 GLY C
   residue 1161 GLY C
   residue 1162 PRO C
   pdb C_1163_to_1202.pdb
}
coordpdb C_912_to_988.pdb C
coordpdb C_1163_to_1202.pdb C
coord C 1156 N [cacoIn_nOut 988 C 0]
##### output of python3 parse_pdb_psfgen.py 6m0j.pdb above #####

guesscoord

regenerate angles dihedrals

writepsf "my_${BASEPDB}.psf"
writepdb "unrelaxed.pdb"


mol delete top
mol new my_${BASEPDB}.psf
mol addfile unrelaxed.pdb
set a [atomselect top all]
set molid [molinfo top get id]

set nc 1000
set rcut 3.0
set temperature 3.0
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]
set loopindex 0
set nloops [llength $loops]
foreach l $loops {
  set chain [lindex $l 0]
  puts "Relaxing loop $loopindex ($l) out of $nloops"
  set residueList [[atomselect ${molid} "chain $chain and resid [lindex $l 1] to [lindex $l 2] and name CA"] get residue]
  do_loop_mc ${residueList} ${chain} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] $logid
  set loopindex [expr $loopindex + 1]
}
$a writepdb "my_${BASEPDB}_mcOut.pdb"

mol delete top
mol new my_${BASEPDB}.psf
mol addfile my_${BASEPDB}_mcOut.pdb
set molid [molinfo top get id]
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
$a writepdb "my_${BASEPDB}.pdb"

# clean up
foreach f $LOCALFILES { 
  exec rm $f
}

quit

