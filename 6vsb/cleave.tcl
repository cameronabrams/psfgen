# cleave.tcl
# Cleaves each of chains A, B, and C at the furin site
# A -> A+D
# B -> B+E
# C -> C+F
# cameron f abrams cfa22@drexel.edu
package require psfgen
mol new my_6vsb.psf
mol addfile my_6vsb_mcOut.pdb
foreach c { A B C } d { D E F } {
   [atomselect top "chain $c and resid < 686"] writepdb "${c}_cleaved.pdb"
   set res [atomselect top "chain $c and resid > 685"]
   $res set chain $d
   $res writepdb "${d}_cleaved.pdb"
}
resetpsf
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str

foreach c { A B C } d { D E F } {
   segment $c {
       pdb ${c}_cleaved.pdb
       mutate 682 ARG
       mutate 683 ARG
       mutate 685 ARG
   }
   coordpdb ${c}_cleaved.pdb $c
   segment $d {
       pdb ${d}_cleaved.pdb
   }
   coordpdb ${d}_cleaved.pdb $d
   patch DISU ${a}:131 ${a}:166
   patch DISU ${a}:291 ${a}:301
   patch DISU ${a}:336 ${a}:361
   patch DISU ${a}:379 ${a}:432
   patch DISU ${a}:391 ${a}:525
   patch DISU ${a}:538 ${a}:590
   patch DISU ${a}:617 ${a}:649
   patch DISU ${a}:662 ${a}:671
   patch DISU ${d}:738 ${d}:760
   patch DISU ${d}:743 ${d}:749
   patch DISU ${d}:1032 ${d}:1043
   patch DISU ${d}:1082 ${d}:1126
}
guesscoord
writepsf "my_6vsb.psf"
writepdb "my_6vsb_mcOut.pdb"

