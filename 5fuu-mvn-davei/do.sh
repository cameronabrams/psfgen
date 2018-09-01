#!/bin/bash
echo $CHARMRUN
echo $NAMD2
vmd -dispdev text -e mkpsf.tcl > psfgen1.log
$CHARMRUN +p2 $NAMD2 my_complex_vac_stage0.namd > vac0.log
$CHARMRUN +p2 $NAMD2 my_complex_vac_stage1.namd > vac1.log
echo "Done."

