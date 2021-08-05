#!/bin/bash
# topogromacs control script
# run this AFTER do_test.sh has created the minimized/solvated system
# copies the charmm parameter files to the CWD
$PSFGEN_BASEDIR/scripts/cp_charmm.sh my_6m0j_solv_stage2.namd
echo "Executing topogromacs VMD script..."
vmd -dispdev text -e $PSFGEN_BASEDIR/6m0j/tg.tcl 2>&1 > topogromacs.log
gmx editconf -f tg_my_6m0j_i_needsbox.pdb -o tg_my_6m0j_i.pdb -box `cat tg-cell-nm.in` 2>&1 >> topogromacs.log
echo "Done.  Results in topogromacs.log."
echo "Use for grompp: [ -c tg_my_6m0j_i.pdb -p tg_my_6m0j_i.top ]"
# example grompp for production md
gmx grompp -f nvt.mdp -c tg_my_6m0j_i.pdb -p tg_my_6m0j_i.top -o tg_my_6m0j_i.tpr -maxwarn 2
