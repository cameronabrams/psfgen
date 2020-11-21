if {![info exists PSFGEN_BASEDIR]} {
    if {[info exists env(PSFGEN_BASEDIR)]} {
        set PSFGEN_BASEDIR $env(PSFGEN_BASEDIR)
    } else {
        set PSFGEN_BASEDIR $env(HOME)/research/psfgen
    }
}
set LOCAL_TOPPARDIR $PSFGEN_BASEDIR/charmm
if {![info exists CHARMM_TOPPARDIR]} {
    if {[info exists env(CHARMM_TOPPARDIR)]} {
        set TOPPARDIR $env(CHARMM_TOPPARDIR)
    } else {
        set TOPPARDIR $env(HOME)/charmm/toppar
    }
}
source ${PSFGEN_BASEDIR}/src/loopmc.tcl
source ${PSFGEN_BASEDIR}/scripts/vmdrc.tcl
package require psfgen
topology $TOPPARDIR/top_all36_prot.rtf
topology $TOPPARDIR/stream/carb/toppar_all36_carb_glycopeptide.str
topology $LOCAL_TOPPARDIR/top_all36_carb.rtf
topology $LOCAL_TOPPARDIR/toppar_water_ions.str
topology $LOCAL_TOPPARDIR/mylink.top
#topology patches
#coordpdb segname_A.pdb 

mol new x01_6vxx.psf
mol addfile config2.coor
set m0 [molinfo top get id]

set segnames [lsort -unique [[atomselect $m0 "all"] get segname]]

foreach s $segnames {
  [atomselect top "segname $s"] writepdb SEG$s.pdb
  segment $s {
     pdb SEG$s.pdb
  }
  coordpdb SEG$s.pdb $s
}

source patches.inp
foreach c { A B C } {
  set ilist {185 79 164 688 461 262 488 853 640}
  foreach i $ilist {
    set cc [expr $i - 1]
    set n [expr $i + 1]
    set nn [expr $i + 2]
    patch HEAL $c:$cc $c:$i $c:$n $c:$nn
  } 
}

regenerate angles dihedrals
guesscoord

writepsf "ligated.psf"
writepdb "ligated.pdb"
exit
