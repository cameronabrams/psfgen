package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf

segment A {
  first ACE
  last CT3
  pdb skel.pdb
}

coordpdb skel.pdb A

guesscoord

writepsf "tmp_alad.psf"
writepdb "tmp_alad.pdb"

resetpsf

mol new tmp_alad.pdb
set a [atomselect top all]
$a set resname ALAD
$a set resid 1
$a writepdb tmp_alad_2.pdb

#pdbalias atom ALAD inpdb intop
pdbalias atom ALAD CAY CL
pdbalias atom ALAD HY1 HL1
pdbalias atom ALAD HY2 HL2
pdbalias atom ALAD HY3 HL3
pdbalias atom ALAD CY CLP
pdbalias atom ALAD OY OL
pdbalias atom ALAD N NL
pdbalias atom ALAD HN HL
pdbalias atom ALAD C CRP
pdbalias atom ALAD O OR
pdbalias atom ALAD NT NR
pdbalias atom ALAD HNT HR
pdbalias atom ALAD CAT CR
pdbalias atom ALAD HT1 HR1
pdbalias atom ALAD HT2 HR2
pdbalias atom ALAD HT3 HR3
segment A {
   pdb tmp_alad_2.pdb
}
coordpdb tmp_alad_2.pdb

writepsf my_alad.psf
writepdb my_alad.pdb

exit
