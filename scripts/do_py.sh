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
temperature=310
numsteps=(20000)
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
    cfapdbparse_args=${!i}
  fi
  if [ "${!i}" = "-parser_args_2" ]; then
    i=$((i+1))
    cfapdbparse_args_2=${!i}
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
if [ ${#cfapdbparse_args} -ge 0 ]; then
  $PYTHON3 $PSFGEN_BASEDIR/scripts/cfapdbparse/cfapdbparse.py $cfapdbparse_args -psfgen psfgen1.tcl ${PDB}.pdb
else
  $PYTHON3 $PSFGEN_BASEDIR/scripts/cfapdbparse/cfapdbparse.py -psfgen psfgen1.tcl ${PDB}.pdb 
fi
vmd -dispdev text -e psfgen1.tcl > psfgen1.log
grep writepsf psfgen1.tcl | tail -1 | sed s/writepsf/structure/ > namd_header.1
CURRPSF=`cat namd_header.1|awk '{print $2}'`
grep writepdb psfgen1.tcl | tail -1 | sed s/writepdb/coordinates/ >> namd_header.1

# 3. run NAMD in vacuum
echo "Running namd2 on vacuum system..."
cat namd_header.1 $PSFGEN_BASEDIR/templates/vac.namd | \
    sed s/%OUT%/stage1/g | \
    sed s/%SEED%/${seed}/g | \
    sed s/%TEMPERATURE%/${temperature}/g | > stage1.namd
$CHARMRUN +p1 $NAMD2 stage1.namd > stage1.log
vmdtext -e $PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl -args stage1.coor tmp.pdb
cat tmp_header.pdb tmp.pdb > stage1.pdb
CURRPDB=stage1.pdb
# if second set of cfapdbparse arguments is given...
if [ ${#cfapdbparse_args_2} -ge 0 ]; then
    $PYTHON3 $PSFGEN_BASEDIR/scripts/cfapdbparse/cfapdbparse.py $cfapdbparse_args_2 -psfgen psfgen2.tcl stage1.pdb
    vmd -dispdev text -e psfgen2.tcl > psfgen2.log
    grep writepsf psfgen2.tcl | tail -1 | sed s/writepsf/structure/ > namd_header.2
    grep writepdb psfgen2.tcl | tail -1 | sed s/writepdb/coordinates/ >> namd_header.2
    cat namd_header.2 $PSFGEN_BASEDIR/templates/vac.namd | \
       sed s/%OUT%/stage2/g | \
       sed s/%SEED%/${seed}/g | \
       sed s/%TEMPERATURE%/${temperature}/g | > stage2.namd
    $CHARMRUN +p${NPE} $NAMD2 stage2.namd > stage2.log
    vmd -e $PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl -args stage2.coor tmp.pdb
    cat tmp_header.pdb tmp.pdb > stage2.pdb
    CURRPDB=stage2.pdb
fi

# 4. solvate
# expects my_${PDB}_vac.coor
echo "Generating solvated system..."
vmd -dispdev text -e $PSFGEN_BASEDIR/scripts/solv.tcl -args -psf $CURRPSF -pdb $CURRPDB -outpre solvated  > mysolv.log
CURRPSF=solvated.psf
CURRPDB=solvated.pdb

# 5. run NAMD
echo "structure $CURRPSF" > namd_header.3
echo "coordinates $CURRPDB" >> namd_header.3
cp namd_header.3 namd_header.3-0
firsttimestep=0
ls=`echo "${#numsteps[@]} - 1" | bc`
for s in `seq 0 $ls`; do
    echo "Running namd2 (stage $s) on solvated system..."
    cat namd_header.3-$s $PSFGEN_BASEDIR/templates/solv.namd | \
        sed s/%STAGE%/${s}/g | \
        sed s/%OUT%/solv_stage${s}/g | \
        sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
        sed s/%FIRSTTIMESTEP%/$firsttimestep/g > solv_stage${s}.namd
    $CHARMRUN +p${NPE} $NAMD2 solv_stage${s}.namd > solv_stage${s}.log
    firsttimestep=`echo "100 + $firsttimestep + ${numsteps[$s]}" | bc`
    ss=$i+1
    cp namd_header.3 namd_header.3-$ss
    echo "bincoordinates solv_stage${s}.coor" >> namd_header.3-$ss
    echo "binvelocities solv_stage${s}.vel" >> namd_header.3-$ss
    echo "extendedsystem solv_stage${s}.xsc" >> namd_header.3-$ss
done
firsttimestep=0
cat namd_header.3-$ss $PSFGEN_BASEDIR/templates/solv.namd | \
    sed s/%STAGE%/$ss/g | \
    sed s/%OUT%/solv_prod0/g | \
    sed s/%NUMSTEPS%/10000000/g | \
    sed s/%FIRSTTIMESTEP%/$firsttimestep/g > solv_prod0.namd
 
echo "Done.  Created solv_prod0.namd, solv_stage${s}.coor, solv_stage${s}.vel, and solv_stage${s}.xsc."

