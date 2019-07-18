#!/bin/bash
# master test script for generating an ethanol/water cosolvent system
#
# Launch in a clean directory
#
# change these absolute pathnames to match your system
PDB=gxg
CHARMRUN=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/charmrun
NAMD2=${HOME}/namd/NAMD_2.12_Source/Linux-x86_64-g++/namd2
export PSFGEN_BASEDIR=${HOME}/research/psfgen
ARGC=$#
i=1
PM=0.2 ; # concentration of GXG tripeptide in molar 
WE=0.55 ; # wt fraction of ethanol in water/ethanol mixture
L=50; # box size in angstroms
gxg_pdb=my_gxg_q.pdb
gxg_psf=my_gxg.psf
seed=$RANDOM
FP=GLYP
LP=CTER
densgcc=0.7
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-L" ]; then
    i=$((i+1))
    L=${!i}
  fi
  if [ "${!i}" = "-pm" ]; then
    i=$((i+1))
    PM=${!i}
  fi
  if [ "${!i}" = "-we" ]; then
    i=$((i+1))
    WE=${!i}
  fi
  if [ "${!i}" = "-fp" ]; then
    i=$((i+1))
    FP=${!i}
  fi
  if [ "${!i}" = "-lp" ]; then
    i=$((i+1))
    LP=${!i}
  fi
  if [ "${!i}" = "-gxg_pdb" ]; then
    i=$((i+1))
    gxg_pdb=${!i}
  fi
  if [ "${!i}" = "-gxg_psf" ]; then
    i=$((i+1))
    gxg_psf=${!i}
  fi
  if [ "${!i}" = "-densgcc" ]; then
    i=$((i+1))
    densgcc=${!i}
  fi
  if [ "${!i}" = "-namd2" ]; then
    i=$((i+1))
    NAMD2=${!i}
  fi
  if [ "${!i}" = "-charmrun" ]; then
    i=$((i+1))
    CHARMRUN=${!i}
  fi
  if [ "${!i}" = "-seed" ]; then
    i=$((i+1))
    seed=${!i}
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

z=`cat $gxg_psf | cut -b 10-10,35-44 | grep ^A | awk 'BEGIN{s=0}{s+=$2}END{printf("%.0f\n",s);}'`
echo "Charge: $z"
cp $PSFGEN_BASEDIR/${PDB}/my_eoh_q.pdb .
echo "Generating packmol input files..."
tclsh $PSFGEN_BASEDIR/${PDB}/my_${PDB}_mix_packmolin.tcl -seed $seed -pm $PM -we $WE -L $L -eoh_pdb my_eoh_q.pdb -gxg_pdb $gxg_pdb -z $z -densgcc $densgcc > psfgen_mix_1.log
echo "Running packmol..."
packmol < pm-tmp.in
echo "Making mixture psf..."
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/my_${PDB}_mix_mkpsf.tcl -args -firstpatch $FP -lastpatch $LP > psfgen_mix_2.log

# run NAMD
numsteps=( 100 200 400 800 1600 16900 )
ls=`echo "${#numsteps[@]} - 1" | bc`
firsttimestep=100; # stage-0 minimization
for s in `seq 0 $ls`; do
  echo "Running namd2 (stage $s) on mixture system..."
  cat $PSFGEN_BASEDIR/${PDB}/my_${PDB}_mix_solv_stageN.namd | \
      sed s/%STAGE%/${s}/g | \
      sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
      sed s/%FIRSTTIMESTEP%/$firsttimestep/g > my_${PDB}_mix_solv_stage${s}.namd
  $CHARMRUN +p8 $NAMD2 my_${PDB}_mix_solv_stage${s}.namd > mix_solv_stage${s}.log
  firsttimestep=`echo "$firsttimestep + ${numsteps[$s]}" | bc`
done

grep ^ENERGY: mix_solv_stage?.log | awk '{print $2,$19}' > V.dat
gnuplot $PSFGEN_BASEDIR/${PDB}/V.gp

echo "Done."
