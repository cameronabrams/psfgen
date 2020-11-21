# measure_bonds.tcl -- 
# reads lines from input file, each is CHAIN, RES1, RES2
# measures distance between C on RES1 and CA on RES2
# appends that distance to the lines
# writes all lines back out to file

# first argument is PSF, second is PDB or COOR, third is name of input file

set psf [lindex $argv 0]
set coor [lindex $argv 1]
mol new $psf
mol addfile $coor
set a [atomselect top "all"]
$a set occupancy 0

set fp [open [lindex $argv 2] "r"]
set lines [split [read $fp] \n]
close $fp
set RES1 {}
set RES2 {}
set CH {}
foreach l $lines {
    if { [llength $l] > 0} {}
       lappend CH [lindex $l 0]
       lappend RES1 [lindex $l 1]
       lappend RES2 [lindex $l 2]
    }
}

set bl {}
foreach c $CH i $RES1 j $RES2 {
    set ii [[atomselect top "chain $c and resid $i and name C"] get index]
    set jj [[atomselect top "chain $c and resid $j and name N"] get index]
    [atomselect top "chain $c and resid $j and name N"] set occupancy 1
    lappend bl [measure bond [list $ii $jj]]
}

set fp [open [lindex $argv 2] "w"]
foreach c $CH i $RES1 j $RES2 b $bl {
    puts $fp "$c $i $j [format %.4f $b]"
}
close $fp
$a writepdb "fixed.pdb"
exit
