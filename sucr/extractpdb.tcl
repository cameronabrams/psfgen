mol new my_sucr_i.psf
mol addfile [lindex $argv 0] 


set aglc [atomselect top "resname AGLC"]
$aglc set resid 1
set bfru [atomselect top "resname BFRU"]
$bfru set resid 2

set a [atomselect top "resname AGLC BFRU"]
$a set segname SU
$a set chain S
$a writepdb "my_sucr_su.pdb"

exit
