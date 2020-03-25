# cleave.tcl
# Cleaves each of chains A, B, and C at the furin site
# A -> A+D
# B -> B+E
# C -> C+F
# cameron f abrams cfa22@drexel.edu
package require psfgen
mol new my_6vyb.psf
mol addfile my_6vxx_mcOut.pdb
foreach c { A B C } d { D E F } {
   [atomselect top "protein and chain $c and resid < 686"] writepdb "${c}_protein_cleaved.pdb"
   set res [atomselect top "protein and chain $c and resid > 685"]
   $res set chain $d
   $res writepdb "${d}_protein_cleaved.pdb"
   [atomselect top "segname ${c}S and (resid < 1311 or resid 1321)"] writepdb "${c}S_glycan_cleaved.pdb"
   set res [atomselect top "segname ${c}S and resid > 1310 and not resid 1321"]
   $res set chain $d
   $res writepdb "${d}S_glycan_cleaved.pdb"
}
resetpsf
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str

foreach c { A B C } d { D E F } {
   segment $c {
       pdb ${c}_protein_cleaved.pdb
       mutate 682 ARG
       mutate 683 ARG
       mutate 685 ARG
   }
   coordpdb ${c}_protein_cleaved.pdb $c
   segment $d {
       pdb ${d}_protein_cleaved.pdb
   }
   coordpdb ${d}_protein_cleaved.pdb $d
   segment ${c}S {
       pdb ${c}S_glycan_cleaved.pdb
   }
   coordpdb ${c}S_glycan_cleaved.pdb ${c}S
   segment ${d}S {
       pdb ${d}S_glycan_cleaved.pdb
   }
   coordpdb ${d}S_glycan_cleaved.pdb ${d}S
   patch DISU ${c}:131 ${c}:166
   patch DISU ${c}:291 ${c}:301
   patch DISU ${c}:336 ${c}:361
   patch DISU ${c}:379 ${c}:432
   patch DISU ${c}:391 ${c}:525
   patch DISU ${c}:538 ${c}:590
   patch DISU ${c}:617 ${c}:649
   patch DISU ${c}:662 ${c}:671
   patch DISU ${d}:738 ${d}:760
   patch DISU ${d}:743 ${d}:749
   patch DISU ${d}:1032 ${d}:1043
   patch DISU ${d}:1082 ${d}:1126
}

patch NGLB A:61 AS:1301
patch NGLB A:122 AS:1302
patch NGLB A:234 AS:1303
patch NGLB A:282 AS:1305
patch NGLB A:331 AS:1306
patch NGLB A:343 AS:1307
patch NGLB A:603 AS:1308
patch NGLB A:616 AS:1309
patch NGLB A:657 AS:1310
patch NGLB D:709 DS:1311
patch NGLB D:717 DS:1312
patch NGLB D:801 DS:1314
patch NGLB D:1074 DS:1316
patch NGLB D:1098 DS:1317
patch NGLB D:1134 DS:1319
patch 14bb AS:1303 AS:1304
patch 14bb DS:1312 DS:1313
patch 14bb DS:1314 DS:1315
patch 14bb DS:1317 DS:1318
patch 14bb DS:1319 DS:1320
patch NGLB B:61 BS:1301
patch NGLB B:122 BS:1302
patch NGLB B:165 BS:1319
patch NGLB B:234 BS:1303
patch NGLB B:282 BS:1304
patch NGLB B:331 BS:1305
patch NGLB B:343 BS:1306
patch NGLB B:603 BS:1307
patch NGLB B:616 BS:1308
patch NGLB B:657 BS:1309
patch NGLB E:709 ES:1310
patch NGLB E:717 ES:1311
patch NGLB E:801 ES:1312
patch NGLB E:1074 ES:1314
patch NGLB E:1098 ES:1315
patch NGLB E:1134 ES:1317
patch 14bb ES:1312 ES:1313
patch 14bb ES:1315 ES:1316
patch 14bb ES:1317 ES:1318
patch NGLB C:61 CS:1301
patch NGLB C:122 CS:1302
patch NGLB C:165 CS:1320
patch NGLB C:234 CS:1303
patch NGLB C:282 CS:1304
patch NGLB C:331 CS:1305
patch NGLB C:343 CS:1306
patch NGLB C:603 CS:1307
patch NGLB C:616 CS:1308
patch NGLB C:657 CS:1309
patch NGLB F:709 FS:1310
patch NGLB F:717 FS:1311
patch NGLB F:801 FS:1313
patch NGLB F:1074 FS:1315
patch NGLB F:1098 FS:1316
patch NGLB F:1134 FS:1318
patch 14bb FS:1311 FS:1312
patch 14bb FS:1313 FS:1314
patch 14bb FS:1316 FS:1317
patch 14bb FS:1318 FS:1319

guesscoord
regenerate angles dihedrals

writepsf "my_6vyb.psf"
writepdb "my_6vyb_mcOut.pdb"

