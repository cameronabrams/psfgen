#!/bin/bash
# master test script for generating a solvated system 
#
# single ethanol molecule
#
# Copy this file to a clean directory and launch it.
#
# change these absolute pathnames to match your system
SPDB=3tod
PDB=eoh
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
  if [ "${!i}" = "-spdb" ]; then
    i=$((i+1))
    SPDB=${!i}
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

# 1. download ${SPDB}.pdb if it is not already here
if [ ! -e ${SPDB}.pdb ]; then
  echo "Retrieving ${SPDB}.pdb..."
  wget -q http://www.rcsb.org/pdb/files/${SPDB}.pdb
fi

# 2. make the psf
echo "Generating vacuum system..."
if [ ${#psfgen_args} -ge 0 ]; then
  vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/mkpsf_${PDB}.tcl -args ${psfgen_args[*]} > psfgen1.log
else
  vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/mkpsf_${PDB}.tcl > psfgen1.log
fi

# 3. run NAMD on the single molecule
echo "Running namd2 on vacuum system..."
ln -s $PSFGEN_BASEDIR/${PDB}/my_${PDB}_vac.namd .
$CHARMRUN +p1 $NAMD2 my_${PDB}_vac.namd > vac.log

# 4. solvate
echo "Generating solvated system..."
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv.tcl > psfgen2.log

# 5. run NAMD; staging to avoid patch-grid errors
numsteps=( 100 200 19700 )
ls=`echo "${#numsteps[@]} - 1" | bc`
firsttimestep=100; # stage-0 minimization
for s in `seq 0 $ls`; do
  echo "Running namd2 (stage $s) on solvated system..."
  cat $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv_stageN.namd | \
      sed s/%STAGE%/${s}/g | \
      sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
      sed s/%FIRSTTIMESTEP%/$firsttimestep/g > my_${PDB}_solv_stage${s}.namd
  $CHARMRUN +p8 $NAMD2 my_${PDB}_solv_stage${s}.namd > solv_stage${s}.log
  firsttimestep=`echo "$firsttimestep + ${numsteps[$s]}" | bc`
done

vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/extractpdb.tcl -args sol-stage${s}.coor > extractpdb.log

grep ^ENERGY: solv_stage?.log | awk '{print $2,$19}' > V.dat
gnuplot $PSFGEN_BASEDIR/${PDB}/V.gp

echo "Done."
