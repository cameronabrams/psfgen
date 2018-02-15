#!/bin/bash
# master test script for generating a solvated system 
#
# single uncharged alanine dipeptide
#
# Copy this file to a clean directory and launch it.
#
# change these absolute pathnames to match your system
PDB=alad
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

# 1. generate skeletal pdb
cat > skel.pdb << ENDPDB
ATOM      1  CA  ALA A   1      -7.226  15.907   6.548  1.00  0.00      A    C
ATOM      2  HA  ALA A   1      -7.334  14.857   6.793  1.00  0.00      A    H
ATOM      3  CB  ALA A   1      -6.931  16.009   5.055  1.00  0.00      A    C
ATOM      4  HB1 ALA A   1      -7.778  15.581   4.478  1.00  0.00      A    H
ATOM      5  HB2 ALA A   1      -6.806  17.069   4.744  1.00  0.00      A    H
ATOM      6  HB3 ALA A   1      -6.011  15.449   4.780  1.00  0.00      A    H
END
ENDPDB

# 2. make the psf
echo "Generating vacuum system..."
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/mkpsf_${PDB}.tcl > psfgen1.log

# 3. run NAMD
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

grep ^ENERGY: solv_stage?.log | awk '{print $2,$19}' > V.dat
gnuplot $PSFGEN_BASEDIR/${PDB}/V.gp

echo "Done."
