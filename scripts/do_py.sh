#!/bin/bash
# master test script for generating a solvated system
#
# change these absolute pathnames to match your system
PDB=
VMD=/opt/vmd/1.9.4a38/bin/vmd
CHARMRUN=${HOME}/namd/NAMD_2.13_Source/Linux-x86_64-g++/charmrun
NAMD2=${HOME}/namd/NAMD_2.13_Source/Linux-x86_64-g++/namd2
export PSFGEN_BASEDIR=${HOME}/research/psfgen
export PYTHON3=${HOME}/anaconda3/bin/python3

ARGC=$#
STAGE=0
NPE=8
i=1
seed=$RANDOM
temperature=310
numsteps=(20000)
parser_args=()
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
  if [ "${!i}" = "-solv_stage" ]; then
     i=$((i+1))
     ssl=${!i}
     numsteps=($(echo "$ssl" | tr ',' '\n'))
#     echo "### solvated system will run in ${#numsteps[@]} stages: ${numsteps[@]}"
  fi
  if [ "${!i}" = "-npe" ]; then
    i=$((i+1))
    export NPE=${!i}
  fi
  if [ "${!i}" = "-temperature" ]; then
    i=$((i+1))
    export temperature=${!i}
  fi
  if [ "${!i}" = "-parser_args" ]; then
    i=$((i+1))
    parser_args+=("${!i}")
  fi
  i=$((i+1))
done

nparse=${#parser_args[@]}
echo "Parser will run $nparse times in series"
for t in `seq 0 $((nparse-1))`; do
    echo " -> ${parser_args[$t]}"
done

if [ ${#PDB} -eq 0 ]; then
   echo "ERROR: must provide PDB code using -pdb XXXX"
   exit
fi

# 1. download ${PDB}.pdb if it is not already here
TASK=-1
if [ ! -e ${PDB}.pdb ]; then
  TASK=$((TASK+1))
  echo "TASK $TASK: Retrieving ${PDB}.pdb..."
  wget -q http://www.rcsb.org/pdb/files/${PDB}.pdb
fi
CURRPDB=${PDB}.pdb

# 2. cycles of parsing/relaxing
for pi in `seq 0 $((nparse-1))`; do
   TASK=$((TASK+1))
   CURRPSFGEN=psfgen${TASK}.tcl
   CURRPSFLOG=`echo $CURRPSFGEN | sed s/tcl/log/`
   $PYTHON3 $PSFGEN_BASEDIR/scripts/cfapdbparse/cfapdbparse.py ${parser_args[$pi]} -psfgen ${CURRPSFGEN} ${CURRPDB}
   CURRPSF=`grep ^writepsf ${CURRPSFGEN} | tail -1 | awk '{print $2}'`
   CURRPDB=`grep writepdb ${CURRPSFGEN} | tail -1 | awk '{print $NF}'`
   echo "TASK $TASK: Generating vacuum system ${CURRPSF} + ${CURRPDB}..."
   $VMD -dispdev text -e ${CURRPSFGEN} > ${CURRPSFLOG}
   echo "structure ${CURRPSF}" > namd_header.${TASK}
   echo "coordinates ${CURRPDB}" >> namd_header.${TASK}
   cat namd_header.${TASK} $PSFGEN_BASEDIR/templates/vac.namd | \
       sed s/%OUT%/stage${TASK}/g | \
       sed s/%SEED%/${seed}/g | \
       sed s/%TEMPERATURE%/${temperature}/g > run${TASK}.namd
   echo "        ->  Running namd2 on vacuum system ${CURRPSF}+${CURRPDB}..."
   $CHARMRUN +p${NPE} $NAMD2 run${TASK}.namd > run${TASK}.log
   $VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl -args ${CURRPSF} stage${TASK}.coor tmp.pdb
   cat charmm_header.pdb tmp.pdb > stage${TASK}.pdb
  CURRPDB=stage${TASK}.pdb
done

# solvate
TASK=$((TASK+1))
echo "TASK $TASK: Generating solvated system from ${CURRPSF}+${CURRPDB}..."
$VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/solv.tcl -args -psf $CURRPSF -pdb $CURRPDB -outpre stage${TASK}  > mysolv.log
CURRPSF=stage${TASK}.psf
CURRPDB=stage${TASK}.pdb

# 5. run solvated NAMD
TASK=$((TASK+1))
echo "structure $CURRPSF" > namd_header.${TASK}
echo "coordinates $CURRPDB" >> namd_header.${TASK}
cp namd_header.${TASK} namd_header.${TASK}-0
firsttimestep=0
ls=`echo "${#numsteps[@]} - 1" | bc`
for s in `seq 0 $ls`; do
    echo "          -> Running namd2 (stage $s) on solvated system..."
    cat namd_header.${TASK}-$s $PSFGEN_BASEDIR/templates/solv.namd | \
        sed s/%STAGE%/${s}/g | \
        sed s/%OUT%/run${TASK}_stage${s}/g | \
        sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
        sed s/%SEED%/${seed}/g | \
        sed s/%TEMPERATURE%/${temperature}/g | \
        sed s/%FIRSTTIMESTEP%/$firsttimestep/g > run${TASK}_stage${s}.namd
    $CHARMRUN +p${NPE} $NAMD2 run${TASK}_stage${s}.namd > run${TASK}_stage${s}.log
    firsttimestep=`echo "100 + $firsttimestep + ${numsteps[$s]}" | bc`
    ss=$((s+1))
    cp namd_header.${TASK} namd_header.${TASK}-$ss
    echo "bincoordinates run${TASK}_stage${s}.coor" >> namd_header.${TASK}-$ss
    echo "binvelocities  run${TASK}_stage${s}.vel"  >> namd_header.${TASK}-$ss
    echo "extendedsystem run${TASK}_stage${s}.xsc"  >> namd_header.${TASK}-$ss
done
firsttimestep=0
cat namd_header.${TASK}-$ss $PSFGEN_BASEDIR/templates/solv.namd | \
    sed s/%STAGE%/$ss/g | \
    sed s/%OUT%/prod/g | \
    sed s/%NUMSTEPS%/10000000/g | \
    sed s/%SEED%/${seed}/g | \
    sed s/%TEMPERATURE%/${temperature}/g | \
    sed s/%FIRSTTIMESTEP%/$firsttimestep/g > prod.namd
 
echo "Done.  Created prod.namd, run${TASK}_stage${s}.coor, run${TASK}_stage${s}.vel, and run${TASK}_stage${s}.xsc."

