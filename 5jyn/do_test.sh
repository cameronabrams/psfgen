#!/bin/bash
# master script for generating an MD system
#
PDB=5jyn
NAMD2=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2
if [ -z "$PSFGEN_BASEDIR" ] ; then
  PSFGEN_BASEDIR=${HOME}/research/psfgen
fi
ARGC=$#
COLVARS_INP=my_${PDB}_colvars_op.inp
RESTART=0
i=1
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-restart" ]; then
    RESTART=1
  fi
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
  if [ "${!i}" = "-colvars_inp" ]; then
    i=$((i+1))
    export COLVARS_INP=${!i}
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

if [ $RESTART = 0 ]; then
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

  # 3. run NAMD
  echo "Running namd2 on vacuum system..."
  ln -s $PSFGEN_BASEDIR/${PDB}/my_${PDB}_vac.namd .
  $CHARMRUN +p8 $NAMD2 my_${PDB}_vac.namd > vac.log

  # 4. solvate
  echo "Generating solvated/membrane system..."
  if [ ${#psfgen_args} -ge 0 ]; then
    vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv_1.tcl -args ${psfgen_args[*]} > psfgen2_1.log
  else
    vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv_1.tcl > psfgen2_1.log
  fi  
  echo "...Running packmol; may take a while. Issue 'tail -f packmol.log' to view progress."
  packmol < pm-tmp.in > packmol.log
fi
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv_2.tcl > psfgen2_2.log
# 5. run NAMD; staging to avoid patch-grid errors
numsteps=( 100 200 400 800 1600 3200 19700 )
ls=`echo "${#numsteps[@]} - 1" | bc` 
pzz=10000
surfacetension=`grep cellbasisvector3 cell.inp | awk '{print $3/100.0}'`
firsttimestep=100; # stage-0 minimization
for s in `seq 0 $ls`; do
  echo "Running namd2 (stage $s) on solvated/membrane system..."
  cat $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv_stageN.namd | \
      sed s/%STAGE%/${s}/g | \
      sed s/%PZZ%/$pzz/g | \
      sed s/%SURFACETENSION%/$surfacetension/g | \
      sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
      sed s/%FIRSTTIMESTEP%/$firsttimestep/g > my_${PDB}_solv_stage${s}.namd
  $CHARMRUN +p8 $NAMD2 my_${PDB}_solv_stage${s}.namd > solv_stage${s}.log
  firsttimestep=`echo "$firsttimestep + ${numsteps[$s]}" | bc`
done

echo "Done."
