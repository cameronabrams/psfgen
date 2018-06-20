mol new my_b529_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "resname B529"]
$a set segname X
$a set chain X
$a writepdb "my_b529_x.pdb"

exit
