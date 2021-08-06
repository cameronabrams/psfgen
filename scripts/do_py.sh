#!/bin/bash
#
# driver for cfapdbparse.py
# c 2020, 2021 cameron f abrams cfa22@drexel.edu
#
# This script implements a workflow that generates an
# equilibrated solvated PSF/PDB/COOR/VEL/XSC datafile
# set and a configuration file for NAMD.  In its 
# simplest usage, the user needs only to specify 
# a four-byte PDB code and this script does the rest.
# 
# The most important arguments are '-pyparser-args'
# which you can learn more about in the help for
# cfapdbparse.py.
#
# The flag '-make-gromacs' must be supplied with a pdb and topology file name 
# for output if the user wishes to create Gromacs input using topoGromacs
#

if [[ -z "${PSFGEN_BASEDIR}" ]]; then
    PSFGEN_BASEDIR=${HOME}/research/psfgen
    if [[ ! -d $PSFGEN_BASEDIR ]]; then
        echo "Error: No PSFGEN_BASEDIR found."
        exit -1
    fi
fi
source $PSFGEN_BASEDIR/scripts/utils.sh

if [[ -z "${PYPARSER}" ]]; then
    PYPARSER=${PSFGEN_BASEDIR}/scripts/cfapdbparse/cfapdbparse.py
    if [[ ! -f $PYPARSER ]]; then
        echo "Error: No PYPARSER found."
        exit -1
    fi
fi

check_command vmd VMD 
check_command namd2 NAMD2 
check_command charmrun CHARMRUN
check_command python3 PYTHON3 

RCSB=https://files.rcsb.org/download

ARGC=$#
PDB=()  # array of input PDB file names
STAGE=0  # indicator of staged equilibration
NPE=16  # number of processors to use in MD
i=1
seed=$RANDOM
temperature=310
numsteps=(20000)
pyparser_args=()
pyparser=$PYPARSER
PRODUCTION_STEPS=10000000
DO_TOPOGROMACS=0
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-pdb" ]; then
    i=$((i+1))
    PDB+=("${!i}")
  fi
  if [ "${!i}" = "-namd2" ]; then
    i=$((i+1))
    NAMD2=${!i}
  fi
  if [ "${!i}" = "-charmrun" ]; then
    i=$((i+1))
    CHARMRUN=${!i}
  fi
  if [ "${!i}" = "-solv-stage-steps" ]; then
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
  if [ "${!i}" = "-pyparser-args" ]; then
    i=$((i+1))
    pyparser_args+=("${!i}")
  fi
  if [ "${!i}" = "-pyparser" ]; then
    i=$((i+1))
    pyparser=${!i}
  fi
  if [ "${!i}" = "-python3-path" ]; then
    i=$((i+1))
    PYTHON3=${!i}
  fi
  if [ "${!i}" = "-production-steps" ]; then
    i=$((i+1))
    PRODUCTION_STEPS=${!i}
  fi
  if [ "${!i}" = "-make-gromacs" ]; then
    i=$((i+1))
    DO_TOPOGROMACS=1
    TG_TOP=${!i}
    i=$((i+1))
    TG_PDB=${!i}
  fi
  i=$((i+1))
done

echo "#### do_py.sh: NAMD/Gromacs system builder -- Cameron F Abrams -- cfa22@drexel.edu"
echo "#### Command:"
echo "#    do_py.sh ${@}"

nparse=${#pyparser_args[@]}
# handle case of a single parser run with no arguments
if (( $nparse == 0 )) ; then
    nparse=1
    pyparser_args=("")
fi
echo "#### PyParser $pyparser will run $nparse time$(ess $nparse) in series"
for t in `seq 0 $((nparse-1))`; do
    if [ ! "${pyparser_args[$t]}" = "" ]; then
        echo "####     PyParser arguments for run #$((t+1)): \"${pyparser_args[$t]}\""
    fi
done

npdb=${#PDB[@]}
if [ $npdb -eq 0 ]; then
   echo "ERROR: must provide at least one PDB code using -pdb XXXX"
   exit
fi

echo "#### The following $npdb PDB file$(ess $npdb) $(isare $npdb) used"
BASEPDB=${PDB[0]}.pdb
AUXPDB=()
echo "####     Base: $BASEPDB"
for p in `seq 1 $((npdb-1))`; do
    echo "####     Auxiliary $p: ${PDB[$p]}"
    AUXPDB+=("${PDB[$p]}")
done

# download pdb's if necessary
TASK=0
for p in `seq 0 $((npdb-1))`; do
    pdb=${PDB[$p]}
    if [ ! -e ${pdb}.pdb ]; then
        echo "Retrieving ${pdb}.pdb..."
        wget -q ${RCSB}/${pdb}.pdb
    fi
done
CURRPDB=$BASEPDB

# cycles of parsing/relaxing
for pi in `seq 0 $((nparse-1))`; do
   TASK=$((TASK+1))
   CURRPSFGEN=psfgen${TASK}.tcl
   $PYTHON3 $PYPARSER ${pyparser_args[$pi]} -pe ${NPE} -postscript ps${TASK}.sh -psfgen ${CURRPSFGEN} ${CURRPDB}
   ./ps${TASK}.sh $TASK
   if [ $? -ne 0 ]; then
       echo "Postscript ps${TASK}.sh failed."
       exit 1
    fi
   read CURRPSF CURRPDB CURRCFG < .tmpvar
done

# solvate
TASK=$((TASK+1))
echo "TASK $TASK: Generating solvated system config${TASK}.psf/.pdb from ${CURRPSF}+${CURRPDB}..."
$VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/solv.tcl -args -psf $CURRPSF -pdb $CURRPDB -outpre config${TASK}  > mysolv.log
CURRPSF=config${TASK}.psf
CURRPDB=config${TASK}.pdb

# run solvated NAMD
TASK=$((TASK+1))
if [ ! -f "$CURRCFG" ]; then
    echo "Error: Previous-stage NAMD configuration file $CURRCFG not found."
fi
# get parameter designations from last config file
grep ^parameters $CURRCFG | grep -v water > _par.inp
echo "#### No binary inputs yet -- this run begins using PDB coordinates" > _bin.inp
firsttimestep=0
ls=`echo "${#numsteps[@]} - 1" | bc`
for s in `seq 0 $ls`; do
    echo "Running namd2 (stage $s of $ls) on solvated system..."
    lastnamd=run${TASK}_stage${s}.namd
    lastsys=config${TASK}_stage${s}
    cat $PSFGEN_BASEDIR/templates/solv.namd | \
        sed "/#### SYSTEM CONFIGURATION FILES BEGIN/r _bin.inp" | \
        sed "/#### SYSTEM CONFIGURATION FILES END/i structure $CURRPSF" | \
        sed "/#### SYSTEM CONFIGURATION FILES END/i coordinates $CURRPDB" | \
        sed "/#### PARAMETER FILES BEGIN/r _par.inp" | \
        sed s/%STAGE%/${s}/g | \
        sed s/%OUT%/${lastsys}/g | \
        sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
        sed s/%SEED%/${seed}/g | \
        sed s/%TEMPERATURE%/${temperature}/g | \
        sed s/%FIRSTTIMESTEP%/$firsttimestep/g > $lastnamd 
    $CHARMRUN +p${NPE} $NAMD2 $lastnamd > run${TASK}_stage${s}.log
    if [ $? -ne 0 ]; then
        echo "NAMD fails.  Check log file run${TASK}_stage${s}.log"
    fi
    firsttimestep=`echo "100 + $firsttimestep + ${numsteps[$s]}" | bc`
    ss=$((s+1))
    
    echo "bincoordinates ${lastsys}.coor" > _bin.inp
    echo "binvelocities  ${lastsys}.vel"  >> _bin.inp
    echo "extendedsystem ${lastsys}.xsc"  >> _bin.inp
done
echo "NAMD task$(ess ${#numsteps[@]}) complete.  Packaging for production."

# Prep for production MD
# copy all charmm parameter files to this directory and
# create 'par.inp', which contains local dir file names
${PSFGEN_BASEDIR}/scripts/cp_charmm.sh $lastnamd
firsttimestep=0
cat $PSFGEN_BASEDIR/templates/prod.namd | \
    sed "/#### SYSTEM CONFIGURATION FILES BEGIN/r _bin.inp" | \
    sed "/#### SYSTEM CONFIGURATION FILES END/i structure $CURRPSF" | \
    sed "/#### SYSTEM CONFIGURATION FILES END/i coordinates $CURRPDB" | \
    sed "/#### PARAMETER FILES BEGIN/r par.inp" | \
    sed s/%OUT%/prod/g | \
    sed s/%NUMSTEPS%/$PRODUCTION_STEPS/g | \
    sed s/%SEED%/${seed}/g | \
    sed s/%TEMPERATURE%/${temperature}/g | \
    sed s/%FIRSTTIMESTEP%/$firsttimestep/g > prod.namd
 
tar zvcf namdprod.tgz $CURRPSF \
                  $CURRPDB \
                  config${TASK}_stage${s}.coor \
                  config${TASK}_stage${s}.vel \
                  config${TASK}_stage${s}.xsc \
                  `cat par.inp | awk '{print $2}'` \
                  prod.namd

rm cell.inp _bin.inp _par.inp *restart*
echo "Done.  Created namdprod.tgz."
if (( $DO_TOPOGROMACS == 1 )); then
    echo "Executing $PSFGEN_BASEDIR/scripts/tg.sh  -psf $CURRPSF -i config${TASK}_stage${s} -top $TG_TOP -pdb $TG_PDB"
    $PSFGEN_BASEDIR/scripts/tg.sh  -psf $CURRPSF -i config${TASK}_stage${s} -top $TG_TOP -pdb $TG_PDB ; #-mdp nvt.mdp -tpr test2.tpr
    tar zvcf gmx.tgz $TG_TOP $TG_PDB
    echo "Done. Created gmx.tgz."
fi
echo "#### do_py.sh is done."