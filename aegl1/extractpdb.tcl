mol new my_aegl1_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "resname AEG DLS1 DLS2 or segname T3"]
$a set segname X
$a set chain X
$a writepdb "my_aegl1_x.pdb"

exit
