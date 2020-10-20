# extract frame from dcd and save as namdbin/xsc
package require pbctools
set psf [lindex $argv 0]
set dcd [lindex $argv 1]
set frm [lindex $argv 2]

set pfx "frm${frm}"

mol new $psf
mol addfile $dcd waitfor all

animate write namdbin ${pfx}.coor beg $frm end $frm top
pbc writexst ${pfx}.xsc -first $frm -last $frm

exit


