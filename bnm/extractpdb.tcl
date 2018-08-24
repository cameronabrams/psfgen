mol new my_bnm_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "segname X"]
$a writepdb "my_bnm_x.pdb"

exit
