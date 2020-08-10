#!/bin/bash
# master test script for generating a solvated system
#
# change these absolute pathnames to match your system
PDB=
CHARMRUN=${HOME}/namd/NAMD_2.13_Source/Linux-x86_64-g++/charmrun
NAMD2=${HOME}/namd/NAMD_2.13_Source/Linux-x86_64-g++/namd2
export PSFGEN_BASEDIR=${HOME}/research/psfgen
export PYTHON3=${HOME}/anaconda3/bin/python3

ARGC=$#
STAGE=0
NPE=8
i=1
seed=$RANDOM
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
  if [ "${!i}" = "-stage" ]; then
     export STAGE=1
  fi
  if [ "${!i}" = "-npe" ]; then
    i=$((i+1))
    export NPE=${!i}
  fi
  if [ "${!i}" = "-cfapdbparse_args" ]; then
    i=$((i+1))
    j=0
    while [ $i -le $ARGC ]; do
      cfapdbparse_args[$j]=${!i}
      i=$((i+1))
      j=$((j+1))
    done
  fi
  i=$((i+1))
done

# 1. download ${PDB}.pdb if it is not already here
if [ ! -e ${PDB}.pdb ]; then
  echo "Retrieving ${PDB}.pdb..."
  wget -q http://www.rcsb.org/pdb/files/${PDB}.pdb
fi

# 2. make the psf
echo "Generating vacuum system..."
if [ ${#cfapdbparse_args_args} -ge 0 ]; then
  $PYTHON3 $PSFGEN_BASEDIR/scripts/cfapdbparse/cfapdbparse.py -pdb ${PDB}.pdb $cfapdbparse_args
else
  $PYTHON3 $PSFGEN_BASEDIR/scripts/cfapdbparse/cfapdbparse.py -pdb ${PDB}.pdb
fi
vmd -dispdev text -e mkpsf.tcl > psfgen1.log

# 3. run NAMD in vacuum
echo "Running namd2 on vacuum system..."
cat $PSFGEN_BASEDIR/templates/my_XXXX_vac.namd | sed s/%PDB%/${PDB}/g > my_${PDB}_vac.namd
$CHARMRUN +p1 $NAMD2 my_${PDB}_vac.namd > vac.log

# 4. solvate
echo "Generating solvated system..."
cat $PSFGEN_BASEDIR/templates/my_XXXX_solv.tcl | sed s/%PDB%/${PDB}/g > mysolv.tcl
vmd -dispdev text -e mysolv.tcl > mysolv.log

# 5. run NAMD
if [ $STAGE -eq 0 ] ; then
  echo "Running namd2 on solvated system..."
  cat $PSFGEN_BASEDIR/templates/my_${PDB}_solv.namd | sed s/%PDB%/$PDB/g | sed s/%SEED%/$seed/g > my_${PDB}_solv.namd
  $CHARMRUN +p${NPE} $NAMD2 my_${PDB}_solv.namd > solv.log
else
  SYSNAME=${PDB}
  if [ -f $PSFGEN_BASEDIR/${SYSNAME}/my_${SYSNAME}_solv_stageN.namd ]; then
    numsteps=( 100 200 19700 )
    ls=`echo "${#numsteps[@]} - 1" | bc`
    firsttimestep=100; # stage-0 minimization
    for s in `seq 0 $ls`; do
      echo "Running namd2 (stage $s) on solvated system..."
      cat $PSFGEN_BASEDIR/${SYSNAME}/my_${SYSNAME}_solv_stageN.namd | \
        sed s/%STAGE%/${s}/g | \
        sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
        sed s/%FIRSTTIMESTEP%/$firsttimestep/g > my_${SYSNAME}_solv_stage${s}.namd
      $CHARMRUN +p${NPE} $NAMD2 my_${SYSNAME}_solv_stage${s}.namd > solv_stage${s}.log
      firsttimestep=`echo "$firsttimestep + ${numsteps[$s]}" | bc`
    done
  else
     echo "Error: staging selected by $PSFGEN_BASEDIR/${SYSNAME}/my_${SYSNAME}_solv_stageN.namd is not found."
     exit
  fi
fi
echo "Done."
