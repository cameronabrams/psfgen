# VMD/psfgen script for generating psf/pdb pair for PDB 2mb5
# myoglobin
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
# load some custom TcL procedures to set coordinates correctly
source ${PSFGEN_BASEDIR}/tcl/loopmc.tcl
set LOCALFILES {}

mol new 2mb5.pdb
set molid [molinfo top get id]

set p [atomselect top protein]
set h [atomselect top heme]
set w [atomselect top "resname DOD"]
set c [atomselect top "resname CMO and altloc A"]

# generate the segments with appropriate resname changes
set hsp [atomselect top "protein and resid 12 36 48 81 82 97 113 116 119"]
$hsp set resname "HSP"
set hsd [atomselect top "protein and resid 24 64 93"]
$hsd set resname "HSD"
set bbh [atomselect top "protein and name D"]
$bbh set name "HN"
set nta [atomselect top "protein and resid 1 and name D1 D2 D3"]
$nta set name [list HT1 HT2 HT3]
set cta [atomselect top "protein and resid 153 and name O OXT"]
$cta set name [list OT1 OT2]
$p writepdb "tmp_prot.pdb"; lappend LOCALFILES tmp_prot.pdb
$h set resname "HEME"
$h set resid 1
$h writepdb "tmp_heme.pdb"; lappend LOCALFILES tmp_heme.pdb
$w set resname "TIP3"
$w writepdb "tmp_water.pdb"; lappend LOCALFILES tmp_water.pdb
$c set resname CO
$c writepdb "tmp_co.pdb"; lappend LOCALFILES tmp_co.pdb

mol delete $molid

package require psfgen
# Use the requisite charmm36 (july 2016) topology files
topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/stream/prot/toppar_all36_prot_heme.str
topology $env(HOME)/charmm/toppar/toppar_water_ions_namd.str

segment A {
    pdb tmp_prot.pdb
}
pdbalias atom ILE HG12 HG11
pdbalias atom ILE HG13 HG12
pdbalias atom ILE CD1 CD
pdbalias atom ILE HD11 HD1
pdbalias atom ILE HD12 HD2
pdbalias atom ILE HD13 HD3
pdbalias atom LEU HB2 HB1
pdbalias atom LEU HB3 HB2
pdbalias atom PRO HB2 HB1
pdbalias atom PRO HB3 HB2
pdbalias atom PRO HG2 HG1
pdbalias atom PRO HG3 HG2
pdbalias atom PRO HD2 HD1
pdbalias atom PRO HD3 HD2
pdbalias atom PHE HB2 HB1
pdbalias atom PHE HB3 HB2
pdbalias atom TRP HB2 HB1
pdbalias atom TRP HB3 HB2
pdbalias atom TRP DE1 HE1
pdbalias atom MET HB2 HB1
pdbalias atom MET HB3 HB2
pdbalias atom MET HG2 HG1
pdbalias atom MET HG3 HG2
pdbalias atom GLY HA2 HA1
pdbalias atom GLY HA3 HA2
pdbalias atom SER HB2 HB1
pdbalias atom SER HB3 HB2
pdbalias atom SER DG HG1
pdbalias atom THR DG1 HG1
pdbalias atom ASN HB2 HB1
pdbalias atom ASN HB3 HB2
pdbalias atom ASN DD22 HD21
pdbalias atom ASN DD21 HD22
pdbalias atom GLN HB2 HB1
pdbalias atom GLN HB3 HB2
pdbalias atom GLN HG2 HG1
pdbalias atom GLN HG3 HG2
pdbalias atom GLN DE21 HE21
pdbalias atom GLN DE22 HE22
pdbalias atom TYR HB2 HB1
pdbalias atom TYR HB3 HB2
pdbalias atom TYR DH HH
pdbalias atom HSP HIS HSP
pdbalias atom HSP HB2 HB1
pdbalias atom HSP HB3 HB2
pdbalias atom HSP DD1 HD1
pdbalias atom HSP DE2 HE2
pdbalias atom HSD HIS HSD
pdbalias atom HSD HB2 HB1
pdbalias atom HSD HB3 HB2
pdbalias atom HSD DD1 HD1
pdbalias atom LYS HB2 HB1
pdbalias atom LYS HB3 HB2
pdbalias atom LYS HG2 HG1
pdbalias atom LYS HG3 HG2
pdbalias atom LYS HD2 HD1
pdbalias atom LYS HD3 HD2
pdbalias atom LYS HE2 HE1
pdbalias atom LYS HE3 HE2
pdbalias atom LYS DZ1 HZ1
pdbalias atom LYS DZ2 HZ2
pdbalias atom LYS DZ3 HZ3
pdbalias atom ARG HB2 HB1
pdbalias atom ARG HB3 HB2
pdbalias atom ARG HG2 HG1
pdbalias atom ARG HG3 HG2
pdbalias atom ARG HD2 HD1
pdbalias atom ARG HD3 HD2
pdbalias atom ARG DE HE
pdbalias atom ARG DH11 HH11
pdbalias atom ARG DH12 HH12
pdbalias atom ARG DH21 HH21
pdbalias atom ARG DH22 HH22
pdbalias atom ASP HB2 HB1
pdbalias atom ASP HB3 HB2
pdbalias atom GLU HB2 HB1
pdbalias atom GLU HB3 HB2
pdbalias atom GLU HG2 HG1
pdbalias atom GLU HG3 HG2
coordpdb tmp_prot.pdb A

segment HEME {
    pdb tmp_heme.pdb
}
pdbalias atom HEME HHA HA
pdbalias atom HEME HHB HB
pdbalias atom HEME HHC HC
pdbalias atom HEME HHD HD
pdbalias atom HEME HAA HAA1
pdbalias atom HEME HAAA HAA2
pdbalias atom HEME HBA HBA1
pdbalias atom HEME HBAA HBA2
pdbalias atom HEME HMA HMA1
pdbalias atom HEME HMAA HMA2
pdbalias atom HEME HMAB HMA3
pdbalias atom HEME HMB HMB1
pdbalias atom HEME HMBA HMB2
pdbalias atom HEME HMBB HMB3
pdbalias atom HEME HBB HBB1
pdbalias atom HEME HBBA HBB2
pdbalias atom HEME HMC HMC1
pdbalias atom HEME HMCA HMC2
pdbalias atom HEME HMCB HMC3
pdbalias atom HEME HBC HBC1
pdbalias atom HEME HBCA HBC2
pdbalias atom HEME HMD HMD1
pdbalias atom HEME HMDA HMD2
pdbalias atom HEME HMDB HMD3
pdbalias atom HEME HAD HAD1
pdbalias atom HEME HADA HAD2
pdbalias atom HEME HBD HBD1
pdbalias atom HEME HBDA HBD2
coordpdb tmp_heme.pdb HEME

segment CO {
    pdb tmp_co.pdb
}
coordpdb tmp_co.pdb CO

segment XWAT {
    auto none
    pdb tmp_water.pdb
}
pdbalias atom TIP3 O OH2
pdbalias atom TIP3 D1 H1
pdbalias atom TIP3 D2 H2
coordpdb tmp_water.pdb XWAT

# introduce the His93-Fe bond
patch PHEM A:93 HEME:1

# no guesscoord is required
# regenerate angles dihedrals ; no need for this

writepsf my_2mb5.psf
writepdb my_2mb5_raw.pdb
resetpsf
# the raw pdb file already has all atoms
# lappend LOCALFILES my_2mb5_raw.pdb

mol new my_2mb5.psf
mol addfile my_2mb5_raw.pdb
set a [atomselect top all]
$a set beta 0
set fix [atomselect top "noh"]
$fix set beta 1
$a writepdb my_2mb5_fix.pdb

# clean up
foreach f $LOCALFILES {
  exec rm $f
}

puts "Generated my_2mb5.psf, my_2mb5_raw.pdb, and my_2mb5_fix.pdb."

quit
