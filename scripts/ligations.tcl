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

#### TOPOLOGY FILE LIST STARTS
#### TOPOLOGY FILE LIST ENDS

set psf [lindex $argv 0]
set coor [lindex $argv 1]
set outpsf [lindex $argv 2]
set outpdb [lindex $argv 3]

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

#### LIGATION LIST STARTS
#### LIGATION LIST ENDS

regenerate angles dihedrals
guesscoord

writepsf $outpsf
writepdb $outpdb
exit
