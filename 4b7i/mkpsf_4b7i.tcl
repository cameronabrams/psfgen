# VMD/psfgen script for generating psf/pdb pair for PDB 4b7i
# Human IgG Fc Bearing Hybrid-type Glycans
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

mol new 4b7i.pdb

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

##### output of python3 parse_pdb_psfgen.py 4b7i.pdb below #####a
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
set segs  { { A  237  444 } 
            { B  238  265 }  { B  273  293 }  { B  298  322 }  { B  329  444 } 
            }
set loops {
            { B  266  272 }  { B  294  297 }  { B  323  328 }  }
[atomselect top "chain A and protein and resid 237 to 444"] writepdb "A_237_to_444.pdb"
segment A {
   pdb A_237_to_444.pdb
}
coordpdb A_237_to_444.pdb A
[atomselect top "chain B and protein and resid 238 to 265"] writepdb "B_238_to_265.pdb"
[atomselect top "chain B and protein and resid 273 to 293"] writepdb "B_273_to_293.pdb"
[atomselect top "chain B and protein and resid 298 to 322"] writepdb "B_298_to_322.pdb"
[atomselect top "chain B and protein and resid 329 to 444"] writepdb "B_329_to_444.pdb"
segment B {
   pdb B_238_to_265.pdb
   residue 266 VAL B
   residue 267 SER B
   residue 268 HSE B
   residue 269 GLU B
   residue 270 ASP B
   residue 271 PRO B
   residue 272 GLU B
   pdb B_273_to_293.pdb
   residue 294 GLU B
   residue 295 GLN B
   residue 296 TYR B
   residue 297 ASN B
   pdb B_298_to_322.pdb
   residue 323 VAL B
   residue 324 SER B
   residue 325 ASN B
   residue 326 LYS B
   residue 327 ALA B
   residue 328 LEU B
   pdb B_329_to_444.pdb
}
coordpdb B_238_to_265.pdb B
coordpdb B_273_to_293.pdb B
coordpdb B_298_to_322.pdb B
coordpdb B_329_to_444.pdb B
coord B 266 N [cacoIn_nOut 265 B 0]
coord B 294 N [cacoIn_nOut 293 B 0]
coord B 323 N [cacoIn_nOut 322 B 0]
set myseg [atomselect top "chain A and resid 1445 to 1452"]
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
$myseg writepdb "AS_1445_to_1452.pdb"
segment AS {
   pdb AS_1445_to_1452.pdb
}
coordpdb AS_1445_to_1452.pdb AS
set myseg [atomselect top "chain B and resid 1445 to 1445"]
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
$myseg writepdb "BI_1445_to_1445.pdb"
segment BI {
   pdb BI_1445_to_1445.pdb
}
coordpdb BI_1445_to_1445.pdb BI
set myseg [atomselect top "chain A and resid 2001 to 2147"]
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
$myseg writepdb "AWX_2001_to_2147.pdb"
segment AWX {
   pdb AWX_2001_to_2147.pdb
}
coordpdb AWX_2001_to_2147.pdb AWX
set myseg [atomselect top "chain B and resid 2001 to 2071"]
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
$myseg writepdb "BWX_2001_to_2071.pdb"
segment BWX {
   pdb BWX_2001_to_2071.pdb
}
coordpdb BWX_2001_to_2071.pdb BWX
patch NGLB A:297 AS:1445
patch 16[axeq 1451 0 A C1 1445][axeq 1445 0 A O6 -1] AS:1445 AS:1451
patch 14[axeq 1446 0 A C1 1445][axeq 1445 0 A O4 -1] AS:1445 AS:1446
patch 14[axeq 1447 0 A C1 1446][axeq 1446 0 A O4 -1] AS:1446 AS:1447
patch 16[axeq 1448 0 A C1 1447][axeq 1447 0 A O6 -1] AS:1447 AS:1448
patch 13[axeq 1452 0 A C1 1447][axeq 1447 0 A O3 -1] AS:1447 AS:1452
patch 16[axeq 1449 0 A C1 1448][axeq 1448 0 A O6 -1] AS:1448 AS:1449
patch 13[axeq 1450 0 A C1 1448][axeq 1448 0 A O3 -1] AS:1448 AS:1450
##### output of python3 parse_pdb_psfgen.py 4b7i.pdb above #####

guesscoord

regenerate angles dihedrals

writepsf "my_4b7i.psf"
writepdb "unrelaxed.pdb"

resetpsf
mol new my_4b7i.psf
mol addfile unrelaxed.pdb
set a [atomselect top all]
set molid [molinfo top get id]
set nc 1000
set rcut 3.0
set temperature 2.5
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]
foreach l $loops {
  set chain [lindex $l 0]
  set residueList [[atomselect ${molid} "chain $chain and resid [lindex $l 1] to [lindex $l 2] and name CA"] get residue]
  do_loop_mc ${residueList} ${chain} ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] $logid
}
$a writepdb "my_4b7i_mcOut.pdb"

mol delete top
mol new my_4b7i.psf
mol addfile my_4b7i_mcOut.pdb
set molid [molinfo top get id]
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
$a writepdb "my_4b7i.pdb"

# clean up
foreach f $LOCALFILES { 
  exec rm $f
}

quit

