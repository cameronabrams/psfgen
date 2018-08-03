# VMD/psfgen script for generating psf/pdb pair for PDB 4h8w
# HIV-1 gp120 clade A/E
# antibody and sCD4 in crystal structure are not included
# NAG's are included
#
# cameron f abrams (c) 2017
# drexel university
# chemical and biological engineering

# check for any arguments
set S375H 0
set H61Y 0
set Q105H 0
set V108I 0
set NIK474_476DMR 0
foreach arg $argv {
  if { $arg == "S375H" } {
     set S375H 1
  }
  if { $arg == "H61Y" } {
     set H61Y 1
  }
  if { $arg == "Q105H"} {
     set Q105H 1
  }
  if { $arg == "V108I"} {
     set V108I 1
  }
  if { $arg == "NIK474-476DMR"} {
     set NIK474_476DMR 1
  }
}
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

mol new 4h8w.pdb

# save heavy-atom coordinates for the four continuous protein fragments
set p1 [atomselect top "chain G and protein and resid 44 to 124 198 to 301"]
$p1 writepdb p1.pdb; lappend LOCALFILES p1.pdb;
set p2 [atomselect top "chain G and protein and resid 324 to 396"]
$p2 writepdb p2.pdb; lappend LOCALFILES p2.pdb;
set p3 [atomselect top "chain G and protein and resid 410 to 459"]
$p3 writepdb p3.pdb; lappend LOCALFILES p3.pdb;
set p4 [atomselect top "chain G and protein and resid 463 to 492"]
$p4 writepdb p4.pdb; lappend LOCALFILES p4.pdb

# compute the positions of the N's for the first residue in each gap loop
set n318pos [cacoIn_nOut 301 G 0]
set n403pos [cacoIn_nOut 396 G 0]
set n460pos [cacoIn_nOut 459 G 0]

# save the waters that are assigned to chain G that are close to the protein or ligand
set w [atomselect top "chain G and water and same residue as within 4.0 of (protein and chain G)"]
$w writepdb "wG.pdb"; lappend LOCALFILES wG.pdb

# save the glycans
set g [atomselect top "chain G and resname NAG"]
$g set resname BGNA
$g writepdb "gG.pdb"; lappend LOCALFILES gG.pdb

mol delete top
package require psfgen
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb_namd_cfa.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str
topology $env(HOME)/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias atom BGNA N2 N
pdbalias atom BGNA C7 C
pdbalias atom BGNA O7 O
pdbalias atom BGNA C8 CT

# now we build the protein, including some
# missing residues
segment G {
  pdb p1.pdb
  residue 318 GLY G
  residue 319 GLY G
  residue 320 SER G
  residue 321 GLY G
  residue 322 SER G
  residue 323 GLY G
  pdb p2.pdb
  residue 403 GLY G
  residue 404 ASN G
  residue 405 GLU G
  residue 406 THR G
  residue 407 MET G
  residue 408 LYS G
  residue 409 GLY G
  pdb p3.pdb
  residue 460 ALA G
  residue 461 ASN G
  residue 462 ASN G
  pdb p4.pdb
  if { $S375H == 1 } {
    puts "MKPSF> mutating 375 from S to H"
    mutate 375 HSD
  }
  if { $H61Y == 1 } {
    puts "MKPSF> mutating 61 from H to Y"
    mutate 61 TYR
  }
  if { $Q105H == 1 } {
    puts "MKPSF> mutating 105 from Q to H"
    mutate 105 HSD
  }
  if { $V108I == 1 } { 
    puts "MKPSF> mutating 108 from V to I"
    mutate 108 ILE
  }
  if { $NIK474_476DMR == 1 } {
    puts "MKPSF> mutating 474-476 from NIK to DMR"
    mutate 474 ASP
    mutate 475 MET
    mutate 476 ARG
  }
}

coordpdb p1.pdb G
coordpdb p2.pdb G
coordpdb p3.pdb G
coordpdb p4.pdb G

segment GS {
  pdb gG.pdb
}
coordpdb gG.pdb GS

# now we build the crystal water segments
pdbalias residue HOH TIP3
pdbalias atom TIP3 O OH2
segment WTXG {
  auto none
  pdb wG.pdb
}
coordpdb wG.pdb WTXG

# Now we create the disulfide bonds
patch DISU G:54 G:74
patch DISU G:119 G:205
patch DISU G:218 G:247
patch DISU G:228 G:239
patch DISU G:296 G:331
patch DISU G:378 G:445
patch DISU G:385 G:418
patch DISU G:395 G:410

# And link the NAGs
patch NGLB G:262 GS:504
patch NGLB G:386 GS:509
patch NGLB G:295 GS:507
patch NGLB G:241 GS:503
patch NGLB G:289 GS:506
patch NGLB G:88 GS:501
patch NGLB G:276 GS:505
patch NGLB G:334 GS:508
patch NGLB G:392 GS:510
patch NGLB G:234 GS:502

# add hydrogens
guesscoord
regenerate angles dihedrals

writepsf "my_4h8w.psf"
writepdb "my_4h8w_mcIn.pdb" ; lappend LOCALFILES my_4h8w_mcIn.pdb

resetpsf
mol new my_4h8w.psf
mol addfile my_4h8w_mcIn.pdb
set a [atomselect top all]
set molid [molinfo top get id]

set nc 1000
set rcut 3.0
set temperature 2.5
set k 10.0
set r0 1.5
set bg [atomselect ${molid} "noh"]
set residueList [[atomselect ${molid} "chain G and resid 318 to 323 and name CA"] get residue]
do_loop_mc ${residueList} G ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] -1
set residueList [[atomselect ${molid} "chain G and resid 403 to 409 and name CA"] get residue]
do_loop_mc ${residueList} G ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] -1
set residueList [[atomselect ${molid} "chain G and resid 460 to 462 and name CA"] get residue]
do_loop_mc ${residueList} G ${molid} ${k} ${r0} ${bg} ${rcut} ${nc} ${temperature} [irand_dom 1000 9999] -1

$a writepdb "my_4h8w_mcOut.pdb"

mol delete $molid

resetpsf
mol load psf my_4h8w.psf pdb my_4h8w_mcOut.pdb
set a [atomselect top all]
set fix [atomselect top "protein and not hydrogen and not (resid 301 to 324 403 to 410 459 to 463)"]
$a set beta 0
$fix set beta 1
set fix [atomselect top "not protein and not water and noh"]
$fix set beta 1
set wat [atomselect top "name OH2"]
$wat set beta 1
if { $S375H == 1 } {
  [atomselect top "resid 375"] set beta 0
}
if { $H61Y == 1 } {
  [atomselect top "resid 61"] set beta 0
}
if { $Q105H == 1 } {
  [atomselect top "resid 105"] set beta 0
}
if { $V108I == 1 } { 
  [atomselect top "resid 108"] set beta 0
}
if { $NIK474_476DMR == 1 } {
  [atomselect top "resid 474 475 476"] set beta 0
}

$a writepdb "my_4h8w_fix.pdb"


exit
