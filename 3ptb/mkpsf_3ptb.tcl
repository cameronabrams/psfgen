# VMD/psfgen script for generating psf/pdb pair for PDB 3ptb
# Serine Proteinase with Benzamidine Ligand
#
# gourav shrivastav (c) 2018
# cameron f abrams
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


set LOCALFILES {}

mol new 3ptb.pdb
set molid [molinfo top get id]

set p [atomselect top protein]
set b [atomselect top "resname BEN"]
set w [atomselect top "resname HOH"]


# generate the segments with appropriate resname changes
$p writepdb "tmp_prot.pdb"; lappend LOCALFILES tmp_prot.pdb
$b writepdb "tmp_ben.pdb"; lappend LOCALFILES tmp_ben.pdb

$w set name OH2
$w set resname TIP3
$w set segname WX
$w set chain WX
$w writepdb "tmp_water.pdb"; lappend LOCALFILES tmp_water.pdb

mol delete $molid


package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd_nonbfixes.str
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology ${PSFGEN_BASEDIR}/charmm/benzamidium.str

pdbalias residue HIS HSD

segment A {
    pdb tmp_prot.pdb
}
pdbalias atom ILE CD1 CD
coordpdb tmp_prot.pdb A


segment BEN {
    pdb tmp_ben.pdb
}
pdbalias atom BEN C6 C6
pdbalias atom BEN C5 C5
pdbalias atom BEN C4 C4
pdbalias atom BEN C3 C3
pdbalias atom BEN C2 C2
pdbalias atom BEN C1 C1
pdbalias atom BEN C  C7
pdbalias atom BEN N2 N2
pdbalias atom BEN N1 N1

coordpdb tmp_ben.pdb BEN


segment WX {
    auto none
    pdb tmp_water.pdb
}

coordpdb tmp_water.pdb WX

guesscoord

patch DISU A:22 A:157
patch DISU A:42  A:58 
patch DISU A:136 A:201
patch DISU A:168 A:182
patch DISU A:191 A:220
patch DISU A:232 A:128


writepsf my_3ptb.psf
writepdb my_3ptb_raw.pdb
resetpsf


mol new my_3ptb.psf
mol addfile my_3ptb_raw.pdb
set a [atomselect top all]
$a set beta 0
set fix [atomselect top "noh"]
$fix set beta 1
$a writepdb my_3ptb_fix.pdb

# clean up
foreach f $LOCALFILES {
  exec rm $f
}

puts "Generated my_3ptb.psf, my_3ptb_raw.pdb, and my_3ptb_fix.pdb."

quit
