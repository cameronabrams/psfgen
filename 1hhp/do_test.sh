#!/bin/bash
# master test script for generating a solvated system: 1hhp
#
# Copy this file to a clean directory and launch it.
#
# change these absolute pathnames to match your system
CHARMRUN=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun
NAMD2=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2
export PSFGEN_BASEDIR=${HOME}/research/psfgen

# 1. download 1hhp.pdb if it is not already here
if [ ! -e 1hhp.pdb ]; then
  echo "Retrieving 1hhp.pdb..."
  wget -q http://www.rcsb.org/pdb/files/1hhp.pdb
fi

# 2. make the psf
vmd -dispdev text -e $PSFGEN_BASEDIR/1hhp/mkpsf_1hhp.tcl > psfgen1.log

# 3. run NAMD
echo "Running namd2 on vacuum system..."
ln -s $PSFGEN_BASEDIR/1hhp/my_1hhp_vac.namd .
$CHARMRUN +p1 $NAMD2 my_1hhp_vac.namd > vac.log

# 4. solvate
vmd -dispdev text -e $PSFGEN_BASEDIR/1hhp/my_1hhp_solv.tcl > psfgen2.log

# 5. run NAMD
echo "Running namd2 on solvated system..."
ln -s $PSFGEN_BASEDIR/1hhp/my_1hhp_solv.namd .
ln -s $PSFGEN_BASEDIR/1hhp/my_1hhp_colvars_op.inp .
$CHARMRUN +p8 $NAMD2 my_1hhp_solv.namd > solv.log

echo "Done."

