mol new my_dls1_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "resname DLS1"]
$a set segname X
$a set chain X
$a writepdb "my_dls1_x.pdb"

exit
