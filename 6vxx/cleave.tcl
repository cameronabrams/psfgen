# cleave.tcl
# Cleaves each of chains A, B, and C at the furin site
# A -> A+D
# B -> B+E
# C -> C+F
# cameron f abrams cfa22@drexel.edu
package require psfgen
mol new my_6vxx.psf
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
   
   patch NGLB ${c}:61 ${c}S:1301
   patch NGLB ${c}:122 ${c}S:1302
   patch NGLB ${c}:165 ${c}S:1321
   patch NGLB ${c}:234 ${c}S:1303
   patch NGLB ${c}:282 ${c}S:1305
   patch NGLB ${c}:331 ${c}S:1306
   patch NGLB ${c}:343 ${c}S:1307
   patch NGLB ${c}:603 ${c}S:1308
   patch NGLB ${c}:616 ${c}S:1309
   patch NGLB ${c}:657 ${c}S:1310
   patch NGLB ${d}:709 ${d}S:1311
   patch NGLB ${d}:717 ${d}S:1312
   patch NGLB ${d}:801 ${d}S:1314
   patch NGLB ${d}:1074 ${d}S:1316
   patch NGLB ${d}:1098 ${d}S:1317
   patch NGLB ${d}:1134 ${d}S:1319
   patch 14bb ${c}S:1303 ${c}S:1304
   patch 14bb ${d}S:1312 ${d}S:1313
   patch 14bb ${d}S:1314 ${d}S:1315
   patch 14bb ${d}S:1317 ${d}S:1318
   patch 14bb ${d}S:1319 ${d}S:1320
}

guesscoord
regenerate angles dihedrals

writepsf "my_6vxx.psf"
writepdb "my_6vxx_mcOut.pdb"

