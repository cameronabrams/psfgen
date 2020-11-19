mol new my_dls2_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "resname DLS2"]
$a set segname X
$a set chain X
$a writepdb "my_dls2_x.pdb"
exit

