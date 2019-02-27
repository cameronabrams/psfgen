#!/bin/bash
# master script for generating an MD system of HIV-1 gp41 6HB
# from PDB 3CP1 with option to extend with MPER/TMD
#
PDB=3cp1
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
  echo "...Running packmol, stage 1: lower leaflet..."
  packmol < pm-stage1.in > packmol-stage1.log
  echo "...Running packmol, stage 2: upper leaflet..."
  packmol < pm-stage2.in > packmol-stage2.log
  echo "...Running packmol, stage 3: water and ions..."
  packmol < pm-stage3.in > packmol-stage3.log

  vmd -dispdev text -e $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv_2.tcl > psfgen2_2.log
fi # RESTART == 0

# 5. run NAMD; staging to avoid patch-grid errors
numsteps=( 100 200 400 800 1600 6400 25600 51200 )
ls=`echo "${#numsteps[@]} - 1" | bc` 
firsttimestep=100; # stage-0 minimization
for s in `seq 0 $ls`; do
  echo "Running namd2 (stage $s of $ls) on solvated/membrane system..."
  cat $PSFGEN_BASEDIR/${PDB}/my_${PDB}_solv_stageN.namd | \
      sed s/%STAGE%/${s}/g | \
      sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
      sed s/%FIRSTTIMESTEP%/$firsttimestep/g > my_${PDB}_solv_stage${s}.namd
  $CHARMRUN +p8 $NAMD2 my_${PDB}_solv_stage${s}.namd > solv_stage${s}.log
  firsttimestep=`echo "$firsttimestep + ${numsteps[$s]}" | bc`
done

# 6. make a plot of the volume of the system as a function of time
grep ^ENERGY: solv_stage?.log | awk '{print $2,$19}' > V.dat
Vmax=`cat V.dat | awk 'BEGIN{m=0}{if (m<$2) m=$2}END{print m}'`
Vmin=`cat V.dat | awk 'BEGIN{m=99999999}{if (m>$2) m=$2}END{print m}'`
cat > tmp.gp << EOF
set term pdfcairo enhanced color fontscale 0.7 lw 1.5
set out "V.pdf"
set encoding iso_8859_1
set border 3
set xtics nomirror
set ytics nomirror
set xlabel "time, 10^3 steps (1 step = 2 fs)"
set ylabel "volume, 10^3 \305^3"
ymin=int($Vmin/1000) - 1
ymax=int($Vmax/1000) + 1
set yr [ymin:ymax]
set ytics ymin,2,ymax
p "V.dat" u (\$1/2000):(\$2/1000.) not w l
EOF
gnuplot tmp.gp
rm tmp.gp

echo "Done."
