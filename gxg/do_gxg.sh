#!/bin/bash
# master test script for generating a solvated system 
#
# single zwitterionic G-X-G tripeptide
#
# Launch this in a clean directory!
#
# change these absolute pathnames to match your system
PDB=gxg
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

# 1. generate skeletal pdb of the first GLY residue only
cat > skel.pdb << ENDPDB
ATOM     21  N   GLY A   1      36.024   5.268  -1.139  1.00  1.00           N
ATOM     22  CA  GLY A   1      35.580   5.412   0.234  1.00  1.00           C
ATOM     23  C   GLY A   1      36.417   6.428   0.981  1.00  1.00           C
ATOM     24  O   GLY A   1      35.933   7.105   1.889  1.00  1.00           O
ATOM     25  H   GLY A   1      36.510   4.462  -1.412  1.00  1.00           H
END
ENDPDB

# 2. generate the my_gxg.pdb and my_gxg.psf for the tripeptide
echo "Generating tripeptide molecule..."
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/mkpsf_${PDB}.tcl -args ${psfgen_args[*]} > psfgen1.log

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

# 6. extract single pdb from solvated run
vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/extractpdb.tcl -args sol-stage${s}.coor > extractpdb.log

grep ^ENERGY: solv_stage?.log | awk '{print $2,$19}' > V.dat
gnuplot $PSFGEN_BASEDIR/${PDB}/V.gp

echo "Done."
