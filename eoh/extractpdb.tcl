mol new my_eoh_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "resname ETOH"]
$a set segname E
$a set chain E
$a writepdb "my_eoh_q.pdb"

exit
