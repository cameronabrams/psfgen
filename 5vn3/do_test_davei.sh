#!/bin/bash
# master test script for generating a solvated system: 5vn3 5vn3 5vn3 5vn3 5vn3 
#
# Copy this file to a clean directory and launch it.
#
# change these absolute pathnames to match your system
PDB=5vn3
SPDB=5f4p
SYSNAME=${PDB}
if [ ! $CHARMRUN ]; then
  CHARMRUN=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun
  NAMD2=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2
fi
if [ ! $PSFGEN_BASEDIR ]; then
  PSFGEN_BASEDIR=${HOME}/research/psfgen
fi
ARGC=$#
COLVARS_INP=my_${PDB}_colvars_op.inp
SEED=12345
RESTART=0
i=1
while [ $i -le $ARGC ] ; do
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
  if [ "${!i}" = "-colvars_inp" ]; then
    i=$((i+1))
    export COLVARS_INP=${!i}
  fi
  if [ "${!i}" = "-seed" ]; then
    i=$((i+1))
    export SEED=${!i}
  fi
  if [ "${!i}" = "-restart" ]; then
    i=$((i+1))
    export RESTART=${!i}
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

DAVEI=0
for item in "${psfgen_args[@]}"; do
  if [ "$item" == "-davei" ]; then
    DAVEI=1
  fi
done

if [ "$RESTART" -lt "1" ] ; then
# 1. download ${PDB}.pdb if it is not already here
  if [ ! -e ${PDB}.pdb ]; then
    echo "Retrieving ${PDB}.pdb..."
    wget -q http://www.rcsb.org/pdb/files/${PDB}.pdb
  fi
  if [ ! -e ${SPDB}.pdb ]; then
    echo "Retrieving ${SPDB}.pdb..."
    wget -q http://www.rcsb.org/pdb/files/${SPDB}.pdb
  fi

  # 2. make the psf
  echo "Generating vacuum system..."
  if [ ${#psfgen_args} -ge "0" ]; then
    vmd -dispdev text -e $PSFGEN_BASEDIR/${SYSNAME}/mkpsf_${SYSNAME}.tcl -args -seed $SEED ${psfgen_args[*]} > psfgen1.log
  else
    vmd -dispdev text -e $PSFGEN_BASEDIR/${SYSNAME}/mkpsf_${SYSNAME}.tcl -args -seed $SEED > psfgen1.log
  fi
fi

# 3a. run a single vacuum NAMD
if [ "$DAVEI" -eq "0" ] ; then
  if [ "$RESTART" -lt "2" ] ; then
    echo "Running namd2 on vacuum system (stage 1 of 1)..."
    cp -f $PSFGEN_BASEDIR/${SYSNAME}/my_${SYSNAME}_vac.namd .
    $CHARMRUN +p8 $NAMD2 my_${SYSNAME}_vac.namd > vac.log
  fi
else
# 3b. run three vacuum NAMD stages for placing linker/spacer/Trp3 of DAVEI
  if [ "$RESTART" -lt "2" ] ; then
    echo "Running namd2 on vacuum system (stage 1 of 3)..."
    cp -f $PSFGEN_BASEDIR/${SYSNAME}/my_${SYSNAME}_vac_stage1.namd .
    $CHARMRUN +p8 $NAMD2 my_${SYSNAME}_vac_stage1.namd > vac_stage1.log
  fi
  if [ "$RESTART" -lt "3" ] ; then
    echo "Running namd2 on vacuum system (stage 2 of 3)..."
    cp -f $PSFGEN_BASEDIR/${SYSNAME}/my_${SYSNAME}_vac_stage2.namd .
    $CHARMRUN +p8 $NAMD2 my_${SYSNAME}_vac_stage2.namd > vac_stage2.log
  fi
  if [ "$RESTART" -lt "4" ] ; then
    echo "Running namd2 on vacuum system (stage 3)..."
    cp -f $PSFGEN_BASEDIR/${SYSNAME}/calc_colvar_forces.tcl .
    cp -f $PSFGEN_BASEDIR/${SYSNAME}/my_${SYSNAME}_vac_stage3.namd .
    $CHARMRUN +p8 $NAMD2 my_${SYSNAME}_vac_stage3.namd > vac_stage3.log
  fi
fi

# 4. solvate
if [ "$RESTART" -lt "5" ] ; then
echo "Generating solvated system..."
if [ -f my_${SYSNAME}_vac_stage3.coor ]; then
  cp -f my_${SYSNAME}_vac_stage3.coor my_${SYSNAME}_vac.coor
fi
vmd -dispdev text -e $PSFGEN_BASEDIR/${SYSNAME}/my_${SYSNAME}_solv.tcl > psfgen2.log

# 5. run NAMD; staging to avoid patch-grid errors
numsteps=( 100 200 19700 )
ls=`echo "${#numsteps[@]} - 1" | bc`
firsttimestep=100; # stage-0 minimization
for s in `seq 0 $ls`; do
  echo "Running namd2 (stage $s) on solvated system..."
  cat $PSFGEN_BASEDIR/${SYSNAME}/my_${SYSNAME}_solv_stageN.namd | \
      sed s/%STAGE%/${s}/g | \
      sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
      sed s/%FIRSTTIMESTEP%/$firsttimestep/g > my_${SYSNAME}_solv_stage${s}.namd
  $CHARMRUN +p8 $NAMD2 my_${SYSNAME}_solv_stage${s}.namd > solv_stage${s}.log
  firsttimestep=`echo "$firsttimestep + ${numsteps[$s]}" | bc`
done
fi

if [ "$RESTART" -ge "5" ]; then
  echo "RESTART is $RESTART, but max is 4, so nothing is done."
fi

echo "Done."
