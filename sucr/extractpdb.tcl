mol new my_sucr_i.psf
mol addfile [lindex $argv 0] 

set a [atomselect top "resname ALGC BFRU"]

$a set segname SU
$a set chain S

set aglc [atomselect top "resname AGLC"]
$aglc set resid 1
set bfru [atomselect top "resname BFRU"]
$bfru set resid 2

$a writepdb "my_sucr_su.pdb"

exit
