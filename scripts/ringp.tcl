source /home/cfa/research/psfgen/src/loopmc.tcl
mol new x01_6vxx.psf
mol addfile config.pdb
check_pierced_rings 0 6 1.5
check_pierced_rings 0 5 1.5
exit
