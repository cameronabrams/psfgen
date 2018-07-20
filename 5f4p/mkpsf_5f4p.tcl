# VMD/psfgen script for generating psf/pdb pair for PDB 5f4p
# core monomeric HIV-1 gp120 with BNM-III-170 bound
#
# cameron f abrams (c) 2018
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
# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
set LOCALFILES {}

mol new 5f4p.pdb
set molid [molinfo top get id]

# save heavy-atom coordinates for the four continuous protein fragments
# found by inspection
set p1 [atomselect top "chain A and protein and resid 49 to 300"]
$p1 writepdb p1.pdb; lappend LOCALFILES p1.pdb
set p2 [atomselect top "chain A and protein and resid 325 to 457"]
$p2 writepdb p2.pdb; lappend LOCALFILES p2.pdb
set p3 [atomselect top "chain A and protein and resid 464 to 490"]
$p3 writepdb p3.pdb; lappend LOCALFILES p3.pdb

# compute the vector positions of the N atom of the first residue of each
# loop that will have to be modeled in, based on the positions of
# the CA, C, and O atoms of the last crystallographically represented residue
set n317pos [cacoIn_nOut 300 A 0]
set n458pos [cacoIn_nOut 457 A 0]

# save the waters that are assigned to chain A that are close to the protein or ligand
set w [atomselect top "chain A and water and same residue as within 4.0 of (protein and chain A)"]
$w writepdb wA.pdb
lappend LOCALFILES wA.pdb

# save glycans
set nag [atomselect top "chain A and resname NAG"] 
$nag set resname BGNA
$nag writepdb nagA.pdb
lappend LOCALFILES nagA.pdb

# save ligand
set bnm [atomselect top "chain A and resname 5VG"]
$bnm set resname BNM3
$bnm set segname L
$bnm writepdb bnm3.pdb
lappend LOCALFILES bnm3.pdb

mol delete $molid

# Now we begin to build the composite PSF and PDB

package require psfgen

# Use the requisite charmm36 (july 2016) topology files
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd.str
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology $PSFGEN_BASEDIR/charmm/bnm.str

pdbalias atom ILE CD1 CD

# This alias assumes ALL histidines in the protein 
# are deprotonated at the delta position
pdbalias residue HIS HSD

pdbalias residue NAG BGLC
pdbalias atom BGLC C7 C
pdbalias atom BGLC O7 O
pdbalias atom BGLC C8 CT

# now, we need to fix the atom names in the BNM ligand
foreach nbad [exec grep 5VG 5f4p.pdb | grep -w "5VG A" | grep HETATM | cut -b 13-16] ngood [exec grep ATOM ${PSFGEN_BASEDIR}/charmm/bnm.str | cut -b 6-9 | grep -v ^H] {
  pdbalias atom BNM3 $nbad $ngood
} 

# now we build the protein, including some
# missing residues
segment A {
  pdb p1.pdb
  residue 317 ASN A
  residue 318 GLY A
  residue 319 GLY A
  residue 320 SER A
  residue 321 GLY A
  residue 322 SER A
  residue 323 GLY A
  residue 324 GLY A
  pdb p2.pdb
  residue 458 GLY A
  residue 459 GLY A
  residue 460 ASN A
  residue 461 ASP A
  residue 462 ASP A
  residue 463 ASN A
  pdb p3.pdb
}

segment L {
  pdb bnm3.pdb
}

coordpdb p1.pdb A
coordpdb p2.pdb A
coordpdb p3.pdb A
coordpdb bnm3.pdb L

# properly set the cartesian coordinates
# of the N's of the N-terminal
# missing residues of each missing segment
coord A 317 N $n317pos
coord A 458 N $n458pos

# now we build the glycan segments
segment GA {
  pdb nagA.pdb
}
coordpdb nagA.pdb GA

# now we build the crystal water segments
pdbalias residue HOH TIP3
pdbalias atom TIP3 O OH2
segment WTXA {
  auto none
  pdb wA.pdb
}
coordpdb wA.pdb WTXA

# Now we create the disulfide bonds
patch DISU A:54 A:74
patch DISU A:119 A:205
patch DISU A:218 A:247
patch DISU A:228 A:239
patch DISU A:296 A:331
patch DISU A:378 A:445
patch DISU A:385 A:418

# Add glycan linkages
patch NGLB A:234 GA:502
patch NGLB A:262 GA:503
patch NGLB A:276 GA:504
patch NGLB A:289 GA:505
patch NGLB A:386 GA:506
patch NGLB A:392 GA:507

guesscoord

regenerate angles dihedrals

writepsf "my_5f4p.psf"
writepdb "my_5f4p_rawloops.pdb"
resetpsf

lappend LOCALFILES my_5f4p_rawloops.pdb

mol new my_5f4p.psf
mol addfile my_5f4p_rawloops.pdb
set a [atomselect top all]
set molid [molinfo top get id]

set nc 1000
set rcut 3.0
set temperature 2.5
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]                                                                                                                                   
set residueList [[atomselect ${molid} "chain A and resid 317 to 324 and name CA"] get residue]
do_loop_mc ${residueList} A ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999]
set residueList [[atomselect ${molid} "chain A and resid 458 to 463 and name CA"] get residue]
do_loop_mc ${residueList} A ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999]

$a writepdb "my_5f4p_mcOut.pdb"

set fix [atomselect top "protein and noh and not (resid 317 to 324 458 to 463)"]
$a set beta 0
$fix set beta 1
set wat [atomselect top "name OH2"]
$wat set beta 1
set glyhv [atomselect top "glycan and noh"]
$glyhv set beta 1
$a writepdb "my_5f4p_fix.pdb"

# clean up
foreach f $LOCALFILES {
  exec rm $f
}

puts "Generated my_5f4p.psf, my_5f4p_mc.pdb, and my_5f4p_fix.pdb."

exit
