package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_lipid.rtf
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd.str

readpsf my_5jyn.psf
coordpdb prot.pdb

mol new my_5jyn_packed.pdb ; # output of packmol
set l [atomselect top "lipid and same residue as (name N and z>0)"]
$l set segname U
$l set chain U
$l writepdb U.pdb; lappend LOCALFILES U.pdb
set l [atomselect top "lipid and same residue as (name N and z<0)"]
$l set segname L
$l set chain L
$l writepdb L.pdb; lappend LOCALFILES L.pdb
set w [atomselect top "water and same residue as (name OH2 and z>0)"]
$w set segname WT1
$w set chain WT1
$w writepdb WT1.pdb; lappend LOCALFILES WT1.pdb
set w [atomselect top "water and same residue as (name OH2 and z<0)"]
$w set segname WT2
$w set chain WT2
$w writepdb WT2.pdb; lappend LOCALFILES WT2.pdb
set i [atomselect top "name CLA NA and z>0"]
$i set segname I1
$i set chain I1
$i writepdb I1.pdb; lappend LOCALFILES I1.pdb
set i [atomselect top "name CLA NA and z<0"]
$i set segname I2
$i set chain I2
$i writepdb I2.pdb; lappend LOCALFILES I2.pdb
segment U {
  pdb U.pdb
}
segment L {
  pdb L.pdb
}
segment WT1 {
  auto none
  pdb WT1.pdb
}
segment WT2 {
  auto none
  pdb WT2.pdb
}
segment I1 {
  pdb I1.pdb
}
segment I2 {
  pdb I2.pdb
}

coordpdb U.pdb U
coordpdb L.pdb L
coordpdb WT1.pdb WT1
coordpdb WT2.pdb WT2
coordpdb I1.pdb I1
coordpdb I2.pdb I2

regenerate angles dihedrals

set outputname my_5jyn
writepsf "${outputname}_i.psf"
writepdb "${outputname}_i.pdb"

# generate a special PDB with certain residue alpha-carbons tagged for use
# by a colvars module input file for orientational and positional restraints
mol new ${outputname}_i.psf
mol addfile ${outputname}_i.pdb
set a [atomselect top all]
$a set beta 0.0
set b [atomselect top "protein and name CA"]
$b set beta 1.0
$a writepdb ${outputname}_caB1.pdb

foreach f $LOCALFILES {
  exec rm -f $f
}

quit

