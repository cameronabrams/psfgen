# compute rmsd vs frame number
package require pbctools
set psf [lindex $argv 0]
set dcd [lindex $argv 1]
mol new $psf
mol addfile $dcd waitfor all

set a0 [atomselect top "protein and name CA"]
$a0 frame 0
set a [atomselect top "protein and name CA"]
set fp [open "rmsd.dat" "w"]
for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
    $a frame $i
    puts $fp "$i [format %.4f [measure rmsd $a $a0]]"
}
close $fp
exit


