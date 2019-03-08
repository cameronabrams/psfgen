mol new my_gxg_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "resid 1 2 3"]
$a set segname A
$a set chain A
$a writepdb "my_gxg_q.pdb"

exit
