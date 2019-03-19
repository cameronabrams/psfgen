mol new my_gxg_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "residue 0 1 2"]
$a set segname A
$a set chain A
$a writepdb "my_gxg_q.pdb"

exit
