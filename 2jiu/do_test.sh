#!/bin/bash
# master test script for generating a solvated system: ${PDB}
#
# Copy this file to a clean directory and launch it.
#
# change these absolute pathnames to match your system
PDB=2jiu
PDB_1=2gs6
PDB_2=5uds
CHARMRUN=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun
NAMD2=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2
export PSFGEN_BASEDIR=${HOME}/research/psfgen
ARGC=$#
i=1
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-pdb" ]; then
    i=$((i+1))
    PDB=${!i}
  fi
  if [ "${!i}" = "-namd2" ]; then
    i=$((i+1))
    NAMD2=${!i}
  fi
  if [ "${!i}" = "-charmrun" ]; then
    i=$((i+1))
    CHARMRUN=${!i}
  fi
  if [ "${!i}" = "-psfgen_basedir" ]; then
    i=$((i+1))
    export PSFGEN_BASEDIR=${!i}
  fi
  if [ "${!i}" = "-psfgen_args" ]; then
    i=$((i+1))
    j=0
    while [ $i -le $ARGC ]; do
      psfgen_args[$j]=${!i}
      i=$((i+1))
      j=$((j+1))
    done
  fi
  i=$((i+1))
done

# 1. download necessary PDB's
for p in $PDB $PDB_1 $PDB_2; do
  if [ ! -e ${p}.pdb ]; then
    echo "Retrieving ${p}.pdb..."
    wget -q http://www.rcsb.org/pdb/files/${p}.pdb
  fi
done

# 2. make the psf
echo "Generating vacuum system..."
if [ ${#psfgen_args} -ge 0 ]; then
  vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/mkpsf_${PDB}.tcl -args ${psfgen_args[*]} > psfgen1.log
else
  vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/mkpsf_${PDB}.tcl > psfgen1.log
fi

# 3. run NAMD
echo "Running namd2 on vacuum system..."
ln -s $PSFGEN_BASEDIR/${PDB}/my_${PDB}_vac.namd .
$CHARMRUN +p1 $NAMD2 my_${PDB}_vac.namd > vac.log

# 4. solvate
echo "Generating solvated system..."
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv.tcl > psfgen2.log

# 5. run NAMD
echo "Running namd2 on solvated system..."
ln -s $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv.namd .
ln -s $PSFGEN_BASEDIR/${PDB}/my_${PDB}_colvars_op.inp .
$CHARMRUN +p8 $NAMD2 my_${PDB}_solv.namd > solv.log

echo "Done."
