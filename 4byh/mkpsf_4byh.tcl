# VMD/psfgen script for generating psf/pdb pair for PDB 4byh
# Fc with sialic glycan
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
source ${PSFGEN_BASEDIR}/scripts/vmdrc.tcl
set LOCALFILES {}

set DOMC 0

mol new 4byh.pdb

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
pdbalias atom BGNA N2 N
pdbalias residue SIA ANE5
pdbalias atom ANE5 C10 C
pdbalias atom ANE5 C11 CT
pdbalias atom ANE5 N5 N
pdbalias atom ANE5 O1A O11
pdbalias atom ANE5 O1B O12
pdbalias atom ANE5 O10 O

##### output of python3 parse_pdb_psfgen.py 4byh.pdb below #####a
### autoPDBPARSE BEGIN
set RESDICT(ZN) ZN2
set RESDICT(HOH) TIP3
set RESDICT(CL) CLA
set RESDICT(NAG) BGNA
set RESDICT(MAN) AMAN
set RESDICT(BMA) BMAN
set RESDICT(FUC) AFUC
set RESDICT(GAL) BGAL
set RESDICT(SIA) ANE5AC
set ANAMEDICT(CL) CLA
set segs  { { A  238  444 } 
            { B  236  444 }  }
set loops {
            }
[atomselect top "chain A and protein and resid 238 to 444"] writepdb "A_238_to_444.pdb"
segment A {
   pdb A_238_to_444.pdb
}
coordpdb A_238_to_444.pdb A
[atomselect top "chain B and protein and resid 236 to 444"] writepdb "B_236_to_444.pdb"
segment B {
   pdb B_236_to_444.pdb
}
coordpdb B_236_to_444.pdb B
set myseg [atomselect top "chain A and resid 1445 to 1454"]
set sav_nm [$myseg get resname]
set new_nm [list]
foreach r $sav_nm {
   lappend new_nm $RESDICT($r)
}
$myseg set resname $new_nm
set new_nm [list]
set sav_nm [$myseg get name]
foreach r $sav_nm {
   if { [ info exists ANAMEDICT($r) ] } {
      lappend new_nm $ANAMEDICT($r)
   } else {
      lappend new_nm $r
   }
}
$myseg set name $new_nm
$myseg writepdb "AS_1445_to_1454.pdb"
segment AS {
   pdb AS_1445_to_1454.pdb
}
coordpdb AS_1445_to_1454.pdb AS
set myseg [atomselect top "chain A and resid 1455 to 1455"]
set sav_nm [$myseg get resname]
set new_nm [list]
foreach r $sav_nm {
   lappend new_nm $RESDICT($r)
}
$myseg set resname $new_nm
set new_nm [list]
set sav_nm [$myseg get name]
foreach r $sav_nm {
   if { [ info exists ANAMEDICT($r) ] } {
      lappend new_nm $ANAMEDICT($r)
   } else {
      lappend new_nm $r
   }
}
$myseg set name $new_nm
$myseg writepdb "AI_1455_to_1455.pdb"
segment AI {
   pdb AI_1455_to_1455.pdb
}
coordpdb AI_1455_to_1455.pdb AI
set myseg [atomselect top "chain B and resid 1445 to 1454"]
set sav_nm [$myseg get resname]
set new_nm [list]
foreach r $sav_nm {
   lappend new_nm $RESDICT($r)
}
$myseg set resname $new_nm
set new_nm [list]
set sav_nm [$myseg get name]
foreach r $sav_nm {
   if { [ info exists ANAMEDICT($r) ] } {
      lappend new_nm $ANAMEDICT($r)
   } else {
      lappend new_nm $r
   }
}
$myseg set name $new_nm
$myseg writepdb "BS_1445_to_1454.pdb"
segment BS {
   pdb BS_1445_to_1454.pdb
}
coordpdb BS_1445_to_1454.pdb BS
set myseg [atomselect top "chain B and resid 1455 to 1455"]
set sav_nm [$myseg get resname]
set new_nm [list]
foreach r $sav_nm {
   lappend new_nm $RESDICT($r)
}
$myseg set resname $new_nm
set new_nm [list]
set sav_nm [$myseg get name]
foreach r $sav_nm {
   if { [ info exists ANAMEDICT($r) ] } {
      lappend new_nm $ANAMEDICT($r)
   } else {
      lappend new_nm $r
   }
}
$myseg set name $new_nm
$myseg writepdb "BI_1455_to_1455.pdb"
segment BI {
   pdb BI_1455_to_1455.pdb
}
coordpdb BI_1455_to_1455.pdb BI
set myseg [atomselect top "chain A and resid 2001 to 2062"]
set sav_nm [$myseg get resname]
set new_nm [list]
foreach r $sav_nm {
   lappend new_nm $RESDICT($r)
}
$myseg set resname $new_nm
set new_nm [list]
set sav_nm [$myseg get name]
foreach r $sav_nm {
   if { [ info exists ANAMEDICT($r) ] } {
      lappend new_nm $ANAMEDICT($r)
   } else {
      lappend new_nm $r
   }
}
$myseg set name $new_nm
$myseg set name OH2
$myseg writepdb "AWX_2001_to_2062.pdb"
segment AWX {
   pdb AWX_2001_to_2062.pdb
}
coordpdb AWX_2001_to_2062.pdb AWX
set myseg [atomselect top "chain B and resid 2001 to 2104"]
set sav_nm [$myseg get resname]
set new_nm [list]
foreach r $sav_nm {
   lappend new_nm $RESDICT($r)
}
$myseg set resname $new_nm
set new_nm [list]
set sav_nm [$myseg get name]
foreach r $sav_nm {
   if { [ info exists ANAMEDICT($r) ] } {
      lappend new_nm $ANAMEDICT($r)
   } else {
      lappend new_nm $r
   }
}
$myseg set name $new_nm
$myseg set name OH2
$myseg writepdb "BWX_2001_to_2104.pdb"
segment BWX {
   pdb BWX_2001_to_2104.pdb
}
coordpdb BWX_2001_to_2104.pdb BWX
patch NGLB A:297 AS:1445
patch 14[axeq 1446 0 A C1 1445][axeq 1445 0 A O4 -1] AS:1445 AS:1446
patch 14[axeq 1447 0 A C1 1446][axeq 1446 0 A O4 -1] AS:1446 AS:1447
patch 16[axeq 1450 0 A C1 1447][axeq 1447 0 A O6 -1] AS:1447 AS:1450
patch 13[axeq 1448 0 A C1 1447][axeq 1447 0 A O3 -1] AS:1447 AS:1448
patch 12[axeq 1449 0 A C1 1448][axeq 1448 0 A O2 -1] AS:1448 AS:1449
patch 12[axeq 1451 0 A C1 1450][axeq 1450 0 A O2 -1] AS:1450 AS:1451
patch 14[axeq 1452 0 A C1 1451][axeq 1451 0 A O4 -1] AS:1451 AS:1452
patch SA26E AS:1452 AS:1453
patch 16[axeq 1454 0 A C1 1445][axeq 1445 0 A O6 -1] AS:1445 AS:1454
patch NGLB B:297 BS:1445
patch 14[axeq 1446 0 B C1 1445][axeq 1445 0 B O4 -1] BS:1445 BS:1446
patch 13[axeq 1448 0 B C1 1447][axeq 1447 0 B O3 -1] BS:1447 BS:1448
patch 16[axeq 1451 0 B C1 1447][axeq 1447 0 B O6 -1] BS:1447 BS:1451
patch 14[axeq 1447 0 B C1 1446][axeq 1446 0 B O4 -1] BS:1446 BS:1447
patch 12[axeq 1449 0 B C1 1448][axeq 1448 0 B O2 -1] BS:1448 BS:1449
patch 14[axeq 1450 0 B C1 1449][axeq 1449 0 B O4 -1] BS:1449 BS:1450
patch 12[axeq 1452 0 B C1 1451][axeq 1451 0 B O2 -1] BS:1451 BS:1452
patch 14[axeq 1453 0 B C1 1452][axeq 1452 0 B O4 -1] BS:1452 BS:1453
patch 16[axeq 1454 0 B C1 1445][axeq 1445 0 B O6 -1] BS:1445 BS:1454
##### output of python3 parse_pdb_psfgen.py 4byh.pdb above #####

guesscoord

regenerate angles dihedrals

writepsf "my_4byh.psf"
writepdb "unrelaxed.pdb"

mol delete top
mol new my_4byh.psf
mol addfile unrelaxed.pdb
set molid [molinfo top get id]
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
$a writepdb "my_4byh.pdb"

# clean up
foreach f $LOCALFILES { 
  exec rm $f
}

quit

