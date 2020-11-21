# cleaves chains
# not yet implemented -- do not use!

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

set psf [lindex $argv 0]
set coor [lindex $argv 1]
set outpsf [lindex $argv 2]
set outpdb [lindex $argv 3]
set clv [lrange $argv 4 end]

mol new $psf
mol addfile $coor
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

regenerate angles dihedrals
guesscoord

writepsf $outpsf
writepdb $outpdb
exit
