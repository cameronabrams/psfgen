# VMD/psfgen script for generating psf/pdb pair for a single
# AEG-Linker-Trp3 DAVEI molecule with the AEG aligned onto the
# BMD529 coordinates in PDB 5u7o
#
# cameron f abrams (c) 2018
# cfa22@drexel.edu
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

mol new 5u7o.pdb

set a [atomselect top "resname 83J and name H3 H14 H6 H5 H4 C23 C26 C24 C21 C20 C15 C13 O06 N05 C10 H1 H2 C07 H11 H12 C04 H9 H10 C01 H7 H8 N02 C12 O03 C14 O09 C16 C18 H13 N08 H23 C19 C17 C22 O11 C27 H16 H17 H18 C25 H15 N28 C29"]
$a set chain X
$a set resid 1
$a set resname AEG
$a writepdb "aeg_cryst_frag.pdb"

# now we extract everything up to but not including the azole
# and name them according to the names Avogadro gave the raw AEG

set b529_anlist [list H3 H14 H6 H5 H4 C23 C26 C24 C21 C20 C15 C13 O06 N05 C10 H1 H2 C07 H11 H12 C04 H9 H10 C01 H7 H8 N02 C12 O03 C14 O09 C16 C18 H13 N08 H23 C19 C17 C22 O11 C27 H16 H17 H18 C25 H15 N28 C29]
set aeg_anlist [list H2 H1 H5 H4 H3 C1 C2 C3 C4 C6 C5 C7 O1 N1 C9 H8 H9 C11 H10 H11 C8 H6 H7 C10 H12 H13 N2 C12 O2 C13 O3 C14 C15 H15 N3 H14 C17 C16 C18 O4 C21 H17 H18 H19 C19 H26 N4 C20]
foreach i $b529_anlist j $aeg_anlist {
  puts "$i $j"
  set sel [atomselect top "resname AEG and name $i"]
  $sel set name "x$j"
}
foreach j $aeg_anlist {
  [atomselect top "resname AEG and name x$j"] set name $j
}
[atomselect top "resname AEG and name H2 H1 H5 H4 H3 C1 C2 C3 C4 C6 C5 C7 O1 N1 C9 H8 H9 C11 H10 H11 C8 H6 H7 C10 H12 H13 N2 C12 O2 C13 O3 C14 C15 H15 N3 H14 C17 C16 C18 O4 C21 H17 H18 H19 C19 H26 N4 C20 C22"] writepdb "aeg_frag.pdb"

lappend LOCALFILES aeg_frag.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology $PSFGEN_BASEDIR/charmm/aeg.str
topology $PSFGEN_BASEDIR/charmm/dls1.str
topology $PSFGEN_BASEDIR/charmm/dls2.str
topology $PSFGEN_BASEDIR/charmm/al1p.str

segment X {
   pdb aeg_frag.pdb
}
coordpdb aeg_frag.pdb X

segment Y {
   residue 1 DLS1 Y
}

segment Z {
   residue 1 DLS2 Z
   residue 2 DLS2 Z
   residue 3 DLS2 Z
}

segment T3 {
   first none
   residue 1 ASP T3
   residue 2 LYS T3
   residue 3 TRP T3
   residue 4 ALA T3
   residue 5 SER T3
   residue 6 ILE T3
   residue 7 TRP T3
   residue 8 ASN T3
   residue 9 TRP T3
   last CT2
}

patch AL1P X:1 Y:1
patch LL12 Y:1 Z:1
patch LL22 Z:1 Z:2
patch LL22 Z:2 Z:3
patch LL2P Z:3 T3:1

guesscoord
regenerate angles dihedrals

writepsf "my_aegl1.psf"
writepdb "my_aegl1.pdb"

foreach f $LOCALFILES {
  exec /bin/rm -f $f
}

exit

