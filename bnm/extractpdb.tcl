mol new my_bnm_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "resname BNM3"]
$a set segname X
$a set chain X
$a writepdb "my_bnm_x.pdb"

exit
