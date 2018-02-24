mol new my_mtl_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "resname DMANOL"]
$a set segname Q
$a set chain Q
$a writepdb "my_mtl_q.pdb"

exit
