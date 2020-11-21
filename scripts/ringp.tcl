source /home/cfa/research/psfgen/src/loopmc.tcl
set psf [lindex $argv 0]
set namdbin [lindex $argv 1]
check_pierced_rings 0 6 1.5
check_pierced_rings 0 5 1.5
exit
