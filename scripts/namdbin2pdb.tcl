set psf [lindex $argv 0]
set namdbin [lindex $argv 1]
set pdb [lindex $argv 2]

mol new $psf
mol addfile $namdbin waitfor all
set a [atomselect top all]
$a writepdb $pdb
exit

