# VMD/psfgen script for generating psf/pdb pair for a single 
# BNM-III-170 molecule extracted from the 5f4p PDB entry
# Option to make a DAVEI molecue
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
# Trp3 peptide sequence
set WARHEAD_SEQUENCE [list ASP LYS TRP ALA SER ILE TRP ASN TRP ]
# if DAVEI > 0, this signals to make a davei; value is no. of diethoxyaminoacetate repeats
set DAVEI 0
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-davei" } {
     incr a
     set DAVEI [lindex $argv $a]
  }
} 

source ${PSFGEN_BASEDIR}/src/loopmc.tcl

# start with freshly downloaded PDB
mol new 5f4p.pdb
set a [atomselect top "resname 5VG"]
$a set chain X
$a set resid 1
$a set resname BNM3
$a moveby [vecscale -1 [measure center $a weight mass]]

$a writepdb "bnm.pdb"
lappend LOCALFILES bnm.pdb

package require psfgen

topology $env(HOME)/charmm/toppar/top_all36_prot.rtf
topology $env(HOME)/charmm/toppar/top_all36_na.rtf
topology $env(HOME)/charmm/toppar/top_all36_carb.rtf
topology $env(HOME)/charmm/toppar/top_all36_cgenff.rtf
topology ${PSFGEN_BASEDIR}/charmm/bnm_edited.str
topology ${PSFGEN_BASEDIR}/charmm/dls1.str
topology ${PSFGEN_BASEDIR}/charmm/dls2.str
topology ${PSFGEN_BASEDIR}/charmm/al1p_PEG.str

# correctly alias atom names; names in fresh pdb are bad; names in cgenff-generated bnm.str are good (derived from avogadro)
foreach nbad [exec grep 5VG 5f4p.pdb | grep -w "5VG A" | grep HETATM | cut -b 13-16] ngood [exec grep ATOM ${PSFGEN_BASEDIR}/charmm/bnm_edited.str | cut -b 6-9 | grep -v ^H] {
  pdbalias atom BNM3 $nbad $ngood
} 

# just a placeholder; the whole molecule will be "chain X", but it will have various segments
set x X
# segment "X" is the bnm
segment ${x} {
   pdb bnm.pdb
}

if { [expr $DAVEI > 0] == "1" } {
   # segment X1 is the DLS1 linker/spacer (PEG/azide)
   segment ${x}1 {
     residue 1 DLS1 ${x}
   }
   # segment X2 are the DLS2 aminodiethoxyacetate units
   segment ${x}2 {
     for { set ll 0 } { $ll < $DAVEI } { incr ll } {
        set lll [ expr $ll + 1 ]
        residue $lll DLS2 ${x}
     }
   }
   # segment XT is the warhead peptide (e.g., Trp3) with no N-terminal patch and a neutral C-terminal patch
   segment ${x}T {
     first none
     for { set ll 0 } { $ll <  [llength $WARHEAD_SEQUENCE] } { incr ll } {
       set lll [expr $ll + 1]
       residue $lll [lindex $WARHEAD_SEQUENCE $ll] ${x}
     }
     last CT2
   }
} 
 
coordpdb bnm.pdb X
# bonds in DAVEI
if { [expr $DAVEI > 0] == "1" } {
  patch AL1P ${x}:1 ${x}1:1
  patch LL12 ${x}1:1 ${x}2:1
  for { set nll 1 } { $nll < $DAVEI } { incr nll } {
     patch LL22 ${x}2:${nll} ${x}2:[expr $nll+1]
  }
  patch LL2P ${x}2:${DAVEI} ${x}T:1
}
 
guesscoord
regenerate angles dihedrals
writepsf "my_bnm.psf"
writepdb "my_bnm_tmp.pdb"; lappend LOCALFILES my_bnm_tmp.pdb

resetpsf
mol new my_bnm.psf
mol addfile my_bnm_tmp.pdb waitfor all
if { [ expr $DAVEI > 0 ] == "1" } {
  set sel [atomselect top "segname XT"]
}
fold_alpha_helix top $sel extrabonds-XT.inp
[atomselect top all] writepdb "my_bnm.pdb"

foreach f $LOCALFILES {
  exec /bin/rm -f $f
}


exit

