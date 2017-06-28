#!/bin/bash
# master test script for generating a solvated system: 1hhp
#
# Copy this file to a clean directory and launch it.
#
# change these absolute pathnames to match your system
PDB=2mb5
CHARMRUN=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun
NAMD2=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2
export PSFGEN_BASEDIR=${HOME}/research/psfgen

# 1. download 1hhp.pdb if it is not already here
if [ ! -e ${PDB}.pdb ]; then
  echo "Retrieving ${PDB}.pdb..."
  wget -q http://www.rcsb.org/pdb/files/${PDB}.pdb
fi

# 2. make the psf
echo "Building vacuum system..."
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/mkpsf_${PDB}.tcl > psfgen1.log

# 3. run NAMD
echo "Running namd2 on vacuum system..."
ln -s $PSFGEN_BASEDIR/${PDB}/my_${PDB}_vac.namd .
$CHARMRUN +p1 $NAMD2 my_${PDB}_vac.namd > vac.log

# 4. solvate
echo "Building solvated system..."
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv.tcl > psfgen2.log

# 5. run NAMD
echo "Running namd2 on solvated system..."
ln -s $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv.namd .
ln -s $PSFGEN_BASEDIR/${PDB}/my_${PDB}_colvars_op.inp .
$CHARMRUN +p8 $NAMD2 my_${PDB}_solv.namd > solv.log

echo "Done."

