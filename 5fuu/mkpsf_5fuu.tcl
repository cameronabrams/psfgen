# VMD/psfgen script for generating psf/pdb pair for PDB 5fuu
# soluble, cleaved HIV_(JR-FL) trimer with glycans
#  - two PGT151's (chains H,M,L,N) are deleted
#
# cameron f abrams (c) 2017
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

set MPER_665_to_682 [list LYS TRP ALA SER LEU TRP ASN TRP PHE ASP ILE SER ASN TRP LEU TRP TYR ILE]
set MPER_659_to_682 [list GLN GLU LEU LEU GLU LEU ASP LYS TRP ALA SER LEU TRP ASN TRP PHE ASP ILE SER ASN TRP LEU TRP TYR ILE]
set MPER_EXTEND 0
set TM_683_to_709 [list LYS LEU PHE ILE MET ILE VAL GLY GLY LEU VAL GLY LEU ARG ILE VAL PHE ALA VAL LEU SER ILE VAL ASN ARG VAL ARG ]
set TM_EXTEND 0
set seed 12345
set MAN9 {}
set LOG_DCD 0
set logid -1
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-mper-extend" } {
    set MPER_EXTEND 1
  }
  if { $arg == "-tm-extend" } {
    set TM_EXTEND 1
    if { $MPER_EXTEND == "0" } {
       set MPER_EXTEND 1
    }
  }
  if { $arg == "-seed" } {
    incr a
    set seed [lindex $argv $a]
  }
  if { $arg == "-man9" } { 
    incr a
    lappend MAN9 [lindex $argv $a]
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

set DOMC 1

mol new 5fuu.pdb

# extract contiguous protein segments as individual pdb's
set segs { { A 31 57 } { A 64 137 } { A 150 402 } { A 408 508 } { B 521 664 } 
           { C 31 59 } { C 64 136 } { C 152 401 } { C 409 506 } { D 512 658 }
           { E 31 58 } { E 64 136 } { E 150 404 } { E 410 507 } { F 512 664 } }

set loops { { A 58 63 } { A 138 147 } { A 403 407 } 
            { C 60 63 } { C 137 151 } { C 402 408 }
            { E 59 63 } { E 137 147 } { E 405 409 } }

# compute cartesian coordinates of missing N atoms at the beginning of each 
# loop to be modeled in
set ns {}
set mln { { A 58 57 } { A 138 137 } { A 403 402 } { C 60 59 } { C 137 136 } { C 402 401 } { E 59 58 } { E 137 136 } { E 405 404 } }
foreach m $mln {
  lappend ns [cacoIn_nOut [lindex $m 2] [lindex $m 0] 0]
}

foreach s $segs {
  [atomselect top "protein and chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] writepdb "[lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb"
  lappend LOCALFILES [lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb
}

# save each chain's glycan residues with CHARMM-compatible names
foreach chain { A B C D E F } {
  set g [atomselect top "chain $chain and not protein"] 
  [atomselect top "chain $chain and resname NAG"] set resname BGNA
  [atomselect top "chain $chain and resname MAN"] set resname AMAN
  [atomselect top "chain $chain and resname BMA"] set resname BMAN
  [atomselect top "chain $chain and resname FUC"] set resname AFUC
  [atomselect top "chain $chain and resname GAL"] set resname BGAL
  $g writepdb "${chain}-glycan.pdb"
  lappend LOCALFILES ${chain}-glycan.pdb
}

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

pdbalias residue NAG BGLC
pdbalias atom BGLC C7 C
pdbalias atom BGLC O7 O
pdbalias atom BGLC C8 CT

segment A {
  pdb A_31_to_57.pdb
  residue 58 ALA A
  residue 59 LYS A
  residue 60 ALA A
  residue 61 TYR A
  residue 62 ASP A
  residue 63 THR A
  pdb A_64_to_137.pdb
  residue 138 ASN A
  residue 139 THR A
  residue 140 THR A
  residue 141 ASN A
  residue 142 ASP A
  residue 143 SER A
  residue 144 GLU A
  residue 145 GLY A
  residue 146 THR A
  residue 147 MET A
  pdb A_150_to_402.pdb
  residue 403 GLU A
  residue 404 GLY A
  residue 405 SER A
  residue 406 ASN A
  residue 407 ASN A
  pdb A_408_to_508.pdb
}

segment AS {
  pdb A-glycan.pdb
  # model-build Man9 at N262
  if { [lsearch $MAN9 262] != -1 } {
    residue 1269 AMAN A
    residue 1270 AMAN A
    residue 1271 AMAN A
    residue 1272 AMAN A
  }
  # model-build Man9 at N332
  if { [lsearch $MAN9 332] != -1 } {
    residue 2335 AMAN A
    residue 2336 AMAN A
    residue 2337 AMAN A
    residue 2338 AMAN A
    residue 2339 AMAN A
    residue 2340 AMAN A
    residue 2341 AMAN A
    residue 2342 AMAN A
  }
  # model-building the Man9 at N386
  if { [lsearch $MAN9 386] != -1 } {
    residue 2389 AMAN A
    residue 2390 AMAN A
    residue 2391 AMAN A
    residue 2392 AMAN A
    residue 2393 AMAN A
    residue 2394 AMAN A
    residue 2395 AMAN A
    residue 2396 AMAN A
  }
  # model-building the Man9 at N392
  if { [lsearch $MAN9 392] != -1 } {
    residue 2397 AMAN A
    residue 2398 AMAN A
    residue 2399 AMAN A
    residue 2400 AMAN A
    residue 2401 AMAN A
    residue 2402 AMAN A
    residue 2403 AMAN A
    residue 2404 AMAN A
  }
  # model-build Man9 at N448
  if { [lsearch $MAN9 448] != -1 } {
    residue 2450 BMAN A
    residue 2451 AMAN A
    residue 2452 AMAN A
    residue 2453 AMAN A
    residue 2454 AMAN A 
    residue 2455 AMAN A
    residue 2456 AMAN A
    residue 2457 AMAN A
    residue 2458 AMAN A
  }
}

segment B {
  pdb B_521_to_664.pdb
  if { $MPER_EXTEND == "1" } {
    set lr 0
    for { set r 665 } { $r < 683 } { incr r } {
      residue $r [lindex $MPER_665_to_682 $lr] B
      incr lr
    }
  }
  if { $TM_EXTEND == "1" } {
    set lr 0
    for { set r 683 } { $r < 710 } { incr r } {
      residue $r [lindex $TM_683_to_709 $lr ] B
      incr lr
    }
  }
}

segment BS {
  pdb B-glycan.pdb
}

segment C {
  pdb C_31_to_59.pdb
  residue 60 ALA C
  residue 61 TYR C
  residue 62 ASP C
  residue 63 THR C
  pdb C_64_to_136.pdb
  residue 137 THR C
  residue 138 ASN C
  residue 139 THR C
  residue 140 THR C
  residue 141 ASN C
  residue 142 ASP C
  residue 143 SER C
  residue 144 GLU C
  residue 145 GLY C
  residue 146 THR C
  residue 147 MET C
  residue 150 GLU C ; # incorrectly numbered as 148 in the remark on line 256 of 5fuu.pdb
  residue 151 ARG C ; # incorrectly numbered as 149 in the remark on line 257 of 5fuu.pdb
  pdb C_152_to_401.pdb
  residue 402 THR C
  residue 403 GLU C
  residue 404 GLY C
  residue 405 SER C
  residue 406 ASN C
  residue 407 ASN C
  residue 408 THR C
  pdb C_409_to_506.pdb
}

segment CS {
  pdb C-glycan.pdb
  # model-build Man9 at N332
  if { [lsearch $MAN9 332] != -1 } {
    residue 2335 AMAN C
    residue 2336 AMAN C
    residue 2337 AMAN C
    residue 2338 AMAN C
    residue 2339 AMAN C
    residue 2340 AMAN C
    residue 2341 AMAN C
    residue 2342 AMAN C
  }
  # model-building the Man9 at N262
  if { [lsearch $MAN9 262] != -1 } {
    residue 1269 AMAN C
    residue 1270 AMAN C
    residue 1271 AMAN C
    residue 1272 AMAN C
  }
  # model-building the Man9 at N386
  if { [lsearch $MAN9 386] != -1 } {
    residue 2389 AMAN C
    residue 2390 AMAN C
    residue 2391 AMAN C
    residue 2392 AMAN C
    residue 2393 AMAN C
    residue 2394 AMAN C
    residue 2395 AMAN C
    residue 2396 AMAN C
  }
  # model-building the Man9 at N392
  if { [lsearch $MAN9 392] != -1 } {
    residue 1394 BMAN C
    residue 2397 AMAN C
    residue 2398 AMAN C
    residue 2399 AMAN C
    residue 2400 AMAN C
    residue 2401 AMAN C
    residue 2402 AMAN C
    residue 2403 AMAN C
    residue 2404 AMAN C
  }
  # model-build Man9 at N448
  if { [lsearch $MAN9 448] != -1 } {
    residue 2456 AMAN C
    residue 2457 AMAN C
    residue 2458 AMAN C
  }
}

segment D {
  pdb D_512_to_658.pdb
  if { $MPER_EXTEND == "1" } {
    set lr 0
    for { set r 659 } { $r < 683 } { incr r } {
      residue $r [lindex $MPER_659_to_682 $lr] D
      incr lr
    }
  }
  if { $TM_EXTEND == "1" } {
    set lr 0
    for { set r 683 } { $r < 710 } { incr r } {
      residue $r [lindex $TM_683_to_709 $lr ] D
      incr lr
    }
  }
}

segment DS {
  pdb D-glycan.pdb
}

segment E {
  pdb E_31_to_58.pdb
  residue 59 LYS E
  residue 60 ALA E
  residue 61 TYR E
  residue 62 ASP E
  residue 63 THR E
  pdb E_64_to_136.pdb
  residue 137 THR E
  residue 138 ASN E
  residue 139 THR E
  residue 140 THR E
  residue 141 ASN E
  residue 142 ASP E 
  residue 143 SER E
  residue 144 GLU E
  residue 145 GLY E
  residue 146 THR E
  residue 147 MET E
  pdb E_150_to_404.pdb
  residue 405 SER E
  residue 406 ASN E
  residue 407 ASN E
  residue 408 THR E
  residue 409 GLU E
  pdb E_410_to_507.pdb
}

segment ES {
  pdb E-glycan.pdb
  # model-build Man9 at N332
  if { [lsearch $MAN9 332] != -1 } {
    residue 2335 AMAN E
    residue 2336 AMAN E
    residue 2337 AMAN E
    residue 2338 AMAN E
    residue 2339 AMAN E
    residue 2340 AMAN E
    residue 2341 AMAN E
    residue 2342 AMAN E
  }
  # model-building the Man9 at N262
  if { [lsearch $MAN9 262] != -1 } {
    residue 1269 AMAN E
    residue 1270 AMAN E
    residue 1271 AMAN E
    residue 1272 AMAN E
  }
  # model-building the Man9 at N386
  if { [lsearch $MAN9 386] != -1 } {
    residue 2391 AMAN E
    residue 2392 AMAN E
    residue 2393 AMAN E
    residue 2394 AMAN E
    residue 2395 AMAN E
    residue 2396 AMAN E
    residue 2397 AMAN E
  }
  # model-building the Man9 at N392
  if { [lsearch $MAN9 392] != -1 } {
    residue 2398 AMAN E
    residue 2399 AMAN E
    residue 2400 AMAN E
    residue 2401 AMAN E
    residue 2402 AMAN E
    residue 2403 AMAN E
    residue 2404 AMAN E
    residue 2405 AMAN E
  }
  # model-build Man9 at N448
  if { [lsearch $MAN9 448] != -1 } {
    residue 2456 AMAN E
    residue 2457 AMAN E
    residue 2458 AMAN E
  }
}

segment F {
  pdb F_512_to_664.pdb
  if { $MPER_EXTEND == "1" } {
    set lr 0
    for { set r 665 } { $r < 683 } { incr r } {
      residue $r [lindex $MPER_665_to_682 $lr] F
      incr lr
    }
  }
  if { $TM_EXTEND == "1" } {
    set lr 0
    for { set r 683 } { $r < 710 } { incr r } {
      residue $r [lindex $TM_683_to_709 $lr ] F
      incr lr
    }
  }
}

segment FS {
  pdb F-glycan.pdb
}

foreach s $segs {
  coordpdb [lindex $s 0]_[lindex $s 1]_to_[lindex $s 2].pdb [lindex $s 0]
}

coordpdb A-glycan.pdb AS
coordpdb B-glycan.pdb BS
coordpdb C-glycan.pdb CS
coordpdb D-glycan.pdb DS
coordpdb E-glycan.pdb ES
coordpdb F-glycan.pdb FS

foreach m $mln n $ns {
  puts "coord [lindex $m 0] [lindex $m 1] N $n"
  coord [lindex $m 0] [lindex $m 1] N $n
}

patch DISU A:54 A:74
patch DISU A:119 A:205
patch DISU A:126 A:196
patch DISU A:131 A:157
patch DISU A:218 A:247
patch DISU A:228 A:239
patch DISU A:296 A:331
patch DISU A:378 A:445
patch DISU A:385 A:418
patch DISU B:598 B:604
patch DISU C:54 C:74
patch DISU C:119 C:205
patch DISU C:126 C:196
patch DISU C:131 C:157
patch DISU C:218 C:247
patch DISU C:228 C:239
patch DISU C:296 C:331
patch DISU C:378 C:445
patch DISU C:385 C:418
patch DISU D:598 D:604
patch DISU E:54 E:74
patch DISU E:119 E:205
patch DISU E:126 E:196
patch DISU E:131 E:157
patch DISU E:218 E:247
patch DISU E:228 E:239
patch DISU E:296 E:331
patch DISU E:378 E:445
patch DISU E:385 E:418
patch DISU F:598 F:604

# glycan at N88
patch NGLB A:88 AS:1088
patch 14bb AS:1088 AS:1089
# glycan at N135
patch NGLB A:135 AS:1135
# glycan at N156
patch NGLB A:156 AS:1156
patch 14bb AS:1156 AS:1157
patch 14bb AS:1157 AS:1158
# glycan at N160
patch NGLB A:160 AS:1160
patch 14bb AS:1160 AS:1161
patch 14bb AS:1161 AS:1162
# glycan at N241
patch NGLB A:241 AS:1241
patch 14bb AS:1241 AS:1242
patch 14bb AS:1242 AS:1243
# glycan at N262  Man9
patch NGLB A:262 AS:1262
patch 14bb AS:1262 AS:1263
patch 14bb AS:1263 AS:1264
patch 16ab AS:1264 AS:1265
patch 13ab AS:1264 AS:1267
patch 13ab AS:1265 AS:1266
patch 12aa AS:1267 AS:1268
if { [lsearch $MAN9 262] != -1 } {
  patch 16ab AS:1265 AS:1269
  patch 12aa AS:1269 AS:1270
  patch 12aa AS:1266 AS:1271
  patch 12aa AS:1268 AS:1272
}
# glycan at N276
patch NGLB A:276 AS:1276
patch 14bb AS:1276 AS:1277
patch 14bb AS:1277 AS:1278
# glycan at N295
patch NGLB A:295 AS:1295
patch 14bb AS:1295 AS:1296
# glycan at N301
patch NGLB A:301 AS:1301
patch 14bb AS:1301 AS:1302
# glycan at N332  Man9
patch NGLB A:332 AS:1332
patch 14bb AS:1332 AS:1333
patch 14bb AS:1333 AS:1334
if { [lsearch $MAN9 332] != "-1" } {
  patch 16ab AS:1334 AS:2335
  patch 16ab AS:2335 AS:2338
  patch 12aa AS:2338 AS:2339
  patch 13ab AS:2335 AS:2336
  patch 12aa AS:2336 AS:2337
  patch 13ab AS:1334 AS:2340
  patch 12aa AS:2340 AS:2341
  patch 12aa AS:2341 AS:2342
}
# glycan at N339
patch NGLB A:339 AS:1339
patch 14bb AS:1339 AS:1340
# glycan at N355
patch NGLB A:355 AS:1355
# glycan at N362
patch NGLB A:362 AS:1362
patch 14bb AS:1362 AS:1363
patch 14bb AS:1363 AS:1364
# glycan at N386
patch NGLB A:386 AS:1386
patch 14bb AS:1386 AS:1387
patch 14bb AS:1387 AS:1388
# model-built alpha-mannoses to complete the Man9
if { [lsearch $MAN9 386] != -1 } {
  patch 13ab AS:1388 AS:2389
  patch 12aa AS:2389 AS:2390
  patch 12aa AS:2390 AS:2391
  patch 16ab AS:1388 AS:2392
  patch 13ab AS:2392 AS:2393
  patch 12aa AS:2393 AS:2394
  patch 16ab AS:2392 AS:2395
  patch 12aa AS:2395 AS:2396
}
# glycan at N392
patch NGLB A:392 AS:1392
patch 14bb AS:1392 AS:1393
patch 14bb AS:1393 AS:1394
# model-built alpha-mannoses to complete the Man9
if { [lsearch $MAN9 392] != -1 } {
  patch 13ab AS:1394 AS:2397
  patch 12aa AS:2397 AS:2398
  patch 12aa AS:2398 AS:2399
  patch 16ab AS:1394 AS:2400
  patch 13ab AS:2400 AS:2401
  patch 12aa AS:2401 AS:2402
  patch 16ab AS:2400 AS:2403
  patch 12aa AS:2403 AS:2404
}
# glycan at N397
patch NGLB A:397 AS:1397
# glycan at N448
patch NGLB A:448 AS:1448
patch 14bb AS:1448 AS:1449
if { [lsearch $MAN9 448] != "-1" } {
  patch 14bb AS:1449 AS:2450
  patch 16ab AS:2450 AS:2453
  patch 16ab AS:2453 AS:2455
  patch 12aa AS:2455 AS:2456
  patch 13ab AS:2453 AS:2454
  patch 12aa AS:2454 AS:2457
  patch 13ab AS:2450 AS:2451
  patch 12aa AS:2451 AS:2452
  patch 12aa AS:2452 AS:2458
}
# glycan at 611
patch NGLB B:611 BS:1611
patch 14bb BS:1611 BS:1612
patch 16ab BS:1611 BS:1624
patch 14bb BS:1612 BS:1613
patch 16ab BS:1613 BS:1619
patch 13ab BS:1613 BS:1614
patch 12aa BS:1614 BS:1617
patch 12aa BS:1619 BS:1620
patch 16ab BS:1619 BS:1622
patch 14bb BS:1620 BS:1621
patch 14bb BS:1622 BS:1623
# glycan at N616
patch NGLB B:616 BS:1600
# glycan at N625
patch NGLB B:625 BS:1625
# glycan at N537
patch NGLB B:637 BS:1637
patch 14bb BS:1637 BS:1638
patch 16ab BS:1637 BS:1650
patch 14bb BS:1638 BS:1639
patch 13ab BS:1639 BS:1640
patch 16ab BS:1639 BS:1645
patch 12aa BS:1640 BS:1643
patch 14bb BS:1640 BS:1641
patch 14bb BS:1641 BS:1642
patch 14bb BS:1643 BS:1644
patch 12aa BS:1645 BS:1646
patch 16ab BS:1645 BS:1648
patch 14bb BS:1646 BS:1647

# glycan at N88
patch NGLB C:88 CS:1088
patch 14bb CS:1088 CS:1089
patch 14bb CS:1089 CS:1090
# glycan at N135
patch NGLB C:135 CS:1135
patch 14bb CS:1135 CS:1136
# glycan and N156
patch NGLB C:156 CS:1156
patch 14bb CS:1156 CS:1157
# glycan at N160
patch NGLB C:160 CS:1160
patch 14bb CS:1160 CS:1161
patch 14bb CS:1161 CS:1162
# glycan at N241
patch NGLB C:241 CS:1241
patch 14bb CS:1241 CS:1242
patch 14bb CS:1242 CS:1243
patch 16ab CS:1243 CS:1246
patch 13ab CS:1243 CS:1244
# glycan at N262
patch NGLB C:262 CS:1262
patch 14bb CS:1262 CS:1263
patch 14bb CS:1263 CS:1264
patch 13ab CS:1264 CS:1267
patch 16ab CS:1264 CS:1265
patch 13ab CS:1265 CS:1266
patch 12aa CS:1267 CS:1268
if { [lsearch $MAN9 262] != -1 } {
  patch 16ab CS:1265 CS:1269
  patch 12aa CS:1269 CS:1270
  patch 12aa CS:1266 CS:1271
  patch 12aa CS:1268 CS:1272
}
# glycan at N276
patch NGLB C:276 CS:1276
patch 14bb CS:1276 CS:1277
patch 14bb CS:1277 CS:1278
# glycan at N295
patch NGLB C:295 CS:1295
patch 14bb CS:1295 CS:1296
# glycan at N301
patch NGLB C:301 CS:1301
patch 14bb CS:1301 CS:1302
patch 14bb CS:1302 CS:1303
# glycan at N332
patch NGLB C:332 CS:1332
patch 14bb CS:1332 CS:1333
patch 14bb CS:1333 CS:1334
if { [lsearch $MAN9 332] != "-1" } {
  patch 16ab CS:1334 CS:2335
  patch 16ab CS:2335 CS:2338
  patch 12aa CS:2338 CS:2339
  patch 13ab CS:2335 CS:2336
  patch 12aa CS:2336 CS:2337
  patch 13ab CS:1334 CS:2340
  patch 12aa CS:2340 CS:2341
  patch 12aa CS:2341 CS:2342
}
# glycan at N339
patch NGLB C:339 CS:1339
patch 14bb CS:1339 CS:1340
# glycan at N355
patch NGLB C:355 CS:1355
patch 14bb CS:1355 CS:1356
# glycan at N362
patch NGLB C:362 CS:1362
patch 14bb CS:1362 CS:1363
patch 14bb CS:1363 CS:1364
# glycan at N386
patch NGLB C:386 CS:1386
patch 14bb CS:1386 CS:1387
patch 14bb CS:1387 CS:1388
# model-built alpha-mannoses to complete the Man9
if { [lsearch $MAN9 386] != -1 } {
  patch 13ab CS:1388 CS:2389
  patch 12aa CS:2389 CS:2390
  patch 12aa CS:2390 CS:2391
  patch 16ab CS:1388 CS:2392
  patch 13ab CS:2392 CS:2393
  patch 12aa CS:2393 CS:2394
  patch 16ab CS:2392 CS:2395
  patch 12aa CS:2395 CS:2396
}
# glycan at N392
patch NGLB C:392 CS:1392
patch 14bb CS:1392 CS:1393
# model-built alpha-mannoses to complete the Man9
if { [lsearch $MAN9 392] != -1 } {
  patch 14bb CS:1393 CS:1394
  patch 13ab CS:1394 CS:2397
  patch 12aa CS:2397 CS:2398
  patch 12aa CS:2398 CS:2399
  patch 16ab CS:1394 CS:2400
  patch 13ab CS:2400 CS:2401
  patch 12aa CS:2401 CS:2402
  patch 16ab CS:2400 CS:2403
  patch 12aa CS:2403 CS:2404
}
# glycan at N397
patch NGLB C:397 CS:1397
# glycan at N448
patch NGLB C:448 CS:1448
patch 14bb CS:1448 CS:1449
patch 14bb CS:1449 CS:1450
patch 13ab CS:1450 CS:1451
patch 16ab CS:1450 CS:1453
patch 12aa CS:1451 CS:1452
patch 13ab CS:1453 CS:1454
patch 16ab CS:1453 CS:1455
if { [lsearch $MAN9 448] != "-1" } {
  patch 12aa CS:1455 CS:2456
  patch 12aa CS:1454 CS:2457
  patch 12aa CS:1452 CS:2458
}
# glycan at N611
patch NGLB D:611 DS:1611
patch 14bb DS:1611 DS:1612
patch 16ab DS:1611 DS:1613
# glycan at N616
patch NGLB D:616 DS:1600
# glycan at N625
patch NGLB D:625 DS:1625
patch 14bb DS:1625 DS:1626
# glycan at N637
patch NGLB D:637 DS:1637
patch 14bb DS:1637 DS:1638

# glycan at N88
patch NGLB E:88 ES:1088
patch 14bb ES:1088 ES:1089
patch 14bb ES:1089 ES:1090
# glycan at N135
patch NGLB E:135 ES:1135
patch 14bb ES:1135 ES:1136
# glycan at N156
patch NGLB E:156 ES:1156
patch 14bb ES:1156 ES:1157
patch 14bb ES:1157 ES:1158
# glycan at N160
patch NGLB E:160 ES:1160
patch 14bb ES:1160 ES:1161
# glycan at N187
patch NGLB E:187 ES:1187
# glycan at N241
patch NGLB E:241 ES:1241
patch 14bb ES:1241 ES:1242
patch 14bb ES:1242 ES:1243
patch 16ab ES:1243 ES:1246
patch 13ab ES:1243 ES:1244
patch 12aa ES:1244 ES:1245
# glycan at N262
patch NGLB E:262 ES:1262
patch 14bb ES:1262 ES:1263
patch 14bb ES:1263 ES:1264
patch 13ab ES:1264 ES:1267
patch 16ab ES:1264 ES:1265
patch 13ab ES:1265 ES:1266
patch 12aa ES:1267 ES:1268
if { [lsearch $MAN9 262] != -1 } {
  patch 16ab ES:1265 ES:1269
  patch 12aa ES:1269 ES:1270
  patch 12aa ES:1266 ES:1271
  patch 12aa ES:1268 ES:1272
}
# glycan at N276
patch NGLB E:276 ES:1276
patch 14bb ES:1276 ES:1277
# glycan at N295
patch NGLB E:295 ES:1295
patch 14bb ES:1295 ES:1296
# glycan at N301
patch NGLB E:301 ES:1301
patch 14bb ES:1301 ES:1302
# glycan at N332
patch NGLB E:332 ES:1332
patch 14bb ES:1332 ES:1333
patch 14bb ES:1333 ES:1334
if { [lsearch $MAN9 332] != "-1" } {
  patch 16ab ES:1334 ES:2335
  patch 16ab ES:2335 ES:2338
  patch 12aa ES:2338 ES:2339
  patch 13ab ES:2335 ES:2336
  patch 12aa ES:2336 ES:2337
  patch 13ab ES:1334 ES:2340
  patch 12aa ES:2340 ES:2341
  patch 12aa ES:2341 ES:2342
}
# glycan at N339
patch NGLB E:339 ES:1339
patch 14bb ES:1339 ES:1340
# glycan at N355
patch NGLB E:355 ES:1355
patch 14bb ES:1355 ES:1356
# glycan N362
patch NGLB E:362 ES:1362
patch 14bb ES:1362 ES:1363
# glycan at N386
patch NGLB E:386 ES:1386
patch 14bb ES:1386 ES:1387
patch 14bb ES:1387 ES:1388
patch 16ab ES:1388 ES:1390
# model-built alpha-mannoses to complete the Man9
if { [lsearch $MAN9 386] != -1 } {
  patch 13ab ES:1388 ES:2391
  patch 12aa ES:2391 ES:2392
  patch 12aa ES:2392 ES:2393
  patch 13ab ES:1390 ES:2394
  patch 12aa ES:2394 ES:2395
  patch 16ab ES:1390 ES:2396
  patch 12aa ES:2396 ES:2397
}
# glycan at N392
patch NGLB E:392 ES:1392
patch 14bb ES:1392 ES:1393
patch 14bb ES:1393 ES:1394
# model-built alpha-mannoses to complete the Man9
if { [lsearch $MAN9 392] != -1 } {
  patch 13ab ES:1394 ES:2398
  patch 12aa ES:2398 ES:2399
  patch 12aa ES:2399 ES:2400
  patch 16ab ES:1394 ES:2401
  patch 13ab ES:2401 ES:2402
  patch 12aa ES:2402 ES:2403
  patch 16ab ES:2401 ES:2404
  patch 12aa ES:2404 ES:2405
}
# glycan at N397
patch NGLB E:397 ES:1397
# glycan at N448
patch NGLB E:448 ES:1448
patch 14bb ES:1448 ES:1449
patch 14bb ES:1449 ES:1450
patch 16ab ES:1450 ES:1453
patch 13ab ES:1450 ES:1451
patch 12aa ES:1451 ES:1452
patch 16ab ES:1453 ES:1455
patch 13ab ES:1453 ES:1454
if { [lsearch $MAN9 448] != "-1" } {
  patch 12aa ES:1455 ES:2456
  patch 12aa ES:1454 ES:2457
  patch 12aa ES:1452 ES:2458
}
# glycan at N611
patch NGLB F:611 FS:1611
patch 14bb FS:1611 FS:1612
patch 16ab FS:1611 FS:1624
patch 14bb FS:1612 FS:1613
patch 16ab FS:1613 FS:1619
patch 13ab FS:1613 FS:1614
patch 12aa FS:1614 FS:1617
patch 16ab FS:1619 FS:1622
patch 12aa FS:1619 FS:1620
patch 14bb FS:1620 FS:1621
patch 14bb FS:1622 FS:1623
patch 14bb FS:1625 FS:1626
patch 16ab FS:1637 FS:1650
# glycan at N616
patch NGLB F:616 FS:1600
# glycan at N625
patch NGLB F:625 FS:1625
# glycan at N637
patch NGLB F:637 FS:1637
patch 14bb FS:1637 FS:1638
patch 14bb FS:1638 FS:1639
patch 16ab FS:1639 FS:1645
patch 13ab FS:1639 FS:1640
patch 14bb FS:1640 FS:1641
patch 12aa FS:1640 FS:1643
patch 14bb FS:1641 FS:1642
patch 14bb FS:1643 FS:1644
patch 16ab FS:1645 FS:1648
patch 12aa FS:1645 FS:1646

guesscoord

regenerate angles dihedrals

writepsf "my_5fuu.psf"
writepdb "unrelaxed.pdb"


lappend LOCALFILES unrelaxed.pdb

mol delete top
mol new my_5fuu.psf
mol addfile unrelaxed.pdb
set molid [molinfo top get id]
set or [measure center [atomselect top "all"] weight mass]
set a [atomselect top all]
$a moveby [vecscale -1 $or]
if { $LOG_DCD != "0" } {
   mol new my_5fuu.psf
   mol addfile unrelaxed.pdb
   set logid [molinfo top get id]
   mol top $molid
   log_addframe ${molid} ${logid}
}
set ca [measure center [atomselect top "protein and chain A C E"] weight mass]
set cb [measure center [atomselect top "protein and chain B D F"] weight mass]   
set pi 3.1415928
set dv [vecsub $ca $cb]
set d [veclength $dv]
set cp [expr [lindex $dv 0]/$d]
set sp [expr [lindex $dv 1]/$d]
set p [expr acos($cp)]
if {[expr $sp < 0.0]} {
  set p [expr 2*$pi-$p]
}
set ct [expr [lindex $dv 2]/$d]
set t [expr acos($ct)]
$a move [transaxis z [expr -1 * $p] rad]
$a move [transaxis y [expr -1 * $t] rad]
$a writepdb "unrelaxed2.pdb"
lappend LOCALFILES "unrelaxed2.pdb"

if { $MPER_EXTEND == "1" } {
   foreach b { B D F } {
      set Cterm 682
      if { $TM_EXTEND == "1" } {
         set Cterm 709
      }
      set sel [atomselect $molid "protein and chain $b and resid 658 to $Cterm"]
      fold_alpha_helix $molid $sel 0
      $sel delete
   }
}

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
$a writepdb "my_5fuu_mcOut.pdb"

if { $LOG_DCD != "0" } {
   set loga [atomselect $logid all]
   animate write dcd $log_dcd_file waitfor all sel $loga $logid
}

# make a pdb file that fixes all heavy atoms in the original
# crystal structure -- all added atoms are set as unfixed
# for a minimization
mol delete top
mol new my_5fuu.psf
mol addfile unrelaxed2.pdb
set a [atomselect top all]
$a set beta 0
foreach s $segs {
  [atomselect top "chain [lindex $s 0] and resid [lindex $s 1] to [lindex $s 2]"] set beta 1
}
[atomselect top "not noh"] set beta 0

$a writepdb "my_5fuu_fix.pdb"

# clean up
foreach f $LOCALFILES { 
  exec rm $f
}

quit

