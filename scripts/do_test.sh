#!/bin/bash
# master test script for generating a solvated system: ${PDB}
#
# Copy this file to a clean directory and launch it.
#
# change these absolute pathnames to match your system
PDB=
CHARMRUN=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun
NAMD2=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2
export PSFGEN_BASEDIR=${HOME}/research/psfgen
ARGC=$#
COLVARS_INP=my_${PDB}_colvars_op.inp
STAGE=0
NPE=8
i=1
seed=$RANDOM
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-pdb" ]; then
    i=$((i+1))
    PDB=${!i}
    COLVARS_INP=my_${PDB}_colvars_op.inp
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
  if [ "${!i}" = "-colvars_inp" ]; then
    i=$((i+1))
    export COLVARS_INP=${!i}
  fi
  if [ "${!i}" = "-stage" ]; then
     export STAGE=1
  fi
  if [ "${!i}" = "-npe" ]; then
    i=$((i+1))
    export NPE=${!i}
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

# 1. download ${PDB}.pdb if it is not already here
if [ ! -e ${PDB}.pdb ]; then
  echo "Retrieving ${PDB}.pdb..."
  wget -q http://www.rcsb.org/pdb/files/${PDB}.pdb
fi

# 2. make the psf
echo "Generating vacuum system..."
if [ ${#psfgen_args} -ge 0 ]; then
  vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/mkpsf_${PDB}.tcl -args ${psfgen_args[*]} > psfgen1.log
else
  vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/mkpsf_${PDB}.tcl > psfgen1.log
fi

# 3. run NAMD in vacuum
echo "Running namd2 on vacuum system..."
cp $PSFGEN_BASEDIR/${PDB}/my_${PDB}_vac.namd .
$CHARMRUN +p1 $NAMD2 my_${PDB}_vac.namd > vac.log
if [ -f "go-stage-2" ]; then
   echo "Running namd2 (stage 2, `cat go-stage-2`) on vacuum system..."
   cp $PSFGEN_BASEDIR/${PDB}/my_${PDB}_vac_stage2.namd .
   $CHARMRUN +p1 $NAMD2 my_${PDB}_vac_stage2.namd > vac_stage2.log
else
   for suf in coor vel xsc; do
     if [ -f my_${PDB}_vac_stage1.${suf} ]; then
        mv my_${PDB}_vac_stage1.${suf}  my_${PDB}_vac.${suf}
     fi
   done
fi

# 4. solvate
echo "Generating solvated system..."
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv.tcl > psfgen2.log

# 5. run NAMD
if [ $STAGE -eq 0 ] ; then
  echo "Running namd2 on solvated system..."
  cp $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv.namd .
  cat my_${PDB}_solv.namd | sed s/%SEED%/$seed/g > tmp; mv tmp my_${PDB}_solv.namd
  lcvi=`grep -i colvarsconfig $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv.namd  | awk '{print $2}'`
  if [ -f $PSFGEN_BASEDIR/${PDB}/${COLVARS_INP} ] ; then
     cp $PSFGEN_BASEDIR/${PDB}/${COLVARS_INP} ./${lcvi}
  fi
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
