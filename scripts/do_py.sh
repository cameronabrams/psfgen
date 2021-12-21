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
solvation_log="mysolv.log"
PRODUCTION_STEPS=10000000
DO_TOPOGROMACS=0
SEL_STR=""
CUBICBOX=""
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-pdb" ]; then
    i=$((i+1))
    val="${!i}"
    while [[ $val != -* ]]; do
      PDB+=("$val")
      i=$((i+1))
      if [ $i -gt $ARGC ]; then
         break
      fi
      val="${!i}"
    done
    i=$((i-1))
  elif [ "${!i}" = "-namd2" ]; then
    i=$((i+1))
    NAMD2=${!i}
  elif [ "${!i}" = "-charmrun" ]; then
    i=$((i+1))
    CHARMRUN=${!i}
  elif [ "${!i}" = "-solvation-log" ]; then
    i=$((i+1))
    solvation_log=${!i}
  elif [ "${!i}" = "-solv-stage-steps" ]; then
     i=$((i+1))
     ssl=${!i}
     numsteps=($(echo "$ssl" | tr ',' '\n'))
  elif [ "${!i}" = "-npe" ]; then
    i=$((i+1))
    export NPE=${!i}
  elif [ "${!i}" = "-temperature" ]; then
    i=$((i+1))
    export temperature=${!i}
  elif [ "${!i}" = "-pyparser-args" ]; then
    i=$((i+1))
    pyparser_args+=("${!i}") 
  elif [ "${!i}" = "-pyparser" ]; then
    i=$((i+1))
    pyparser=${!i}
  elif [ "${!i}" = "-python3-path" ]; then
    i=$((i+1))
    PYTHON3=${!i}
  elif [ "${!i}" = "-production-steps" ]; then
    i=$((i+1))
    PRODUCTION_STEPS=${!i}
  elif [ "${!i}" = "-make-gromacs" ]; then
    i=$((i+1))
    DO_TOPOGROMACS=1
    TG_TOP=${!i}
    i=$((i+1))
    TG_PDB=${!i}
  elif [ "${!i}" = "-selection" ]; then
    i=$((i+1))
    SEL_STR=${!i}
  elif [ "${!i}" = "-cubic-box" ]; then
    CUBICBOX="-cubic"
  else
    echo "${!i}: not recognized"
  fi
  i=$((i+1))
done

echo "# do_py.sh: NAMD/Gromacs system builder -- Cameron F Abrams -- cfa22@drexel.edu"
echo "# Uses $pyparser as the PyParser"
echo "# Command:"
echo "#    do_py.sh ${@}"

nparse=${#pyparser_args[@]}
# handle case of a single parser run with no arguments
if (( $nparse == 0 )) ; then
  nparse=1
  pyparser_args=("")
fi
echo "# PyParser $pyparser will run $nparse time$(ess $nparse) in series"
for t in `seq 0 $((nparse-1))`; do
  if [ ! "${pyparser_args[$t]}" = "" ]; then
    echo "#     PyParser arguments for run #$((t+1)): \"${pyparser_args[$t]}\""
  fi
done

npdb=${#PDB[@]}
if [ $npdb -eq 0 ]; then
  echo "ERROR: must provide at least one PDB code using -pdb XXXX"
  exit 1
fi

echo "# The following $npdb PDB file$(ess $npdb) $(isare $npdb) used"
BASEPDB=${PDB[0]}.pdb
AUXPDB=()
echo "#     Base: $BASEPDB"
for p in `seq 1 $((npdb-1))`; do
  echo "#     Auxiliary $p: ${PDB[$p]}"
  AUXPDB+=("${PDB[$p]}")
done

# download pdb's if necessary
TASK=0
for p in `seq 0 $((npdb-1))`; do
  pdb=${PDB[$p]}
  if [ ! -e ${pdb}.pdb ]; then
    echo "# Retrieving ${pdb}.pdb..."
    wget -q ${RCSB}/${pdb}.pdb
  fi
done
CURRPDB=$BASEPDB

# cycles of parsing/relaxing
for pi in `seq 0 $((nparse-1))`; do
  TASK=$((TASK+1))
  CURRPSFGEN=psfgen${TASK}.tcl
  PS=parser-postscript-task${TASK}.sh
  echo "# Task $TASK: Using $pyparser to generate Tcl script and VMD to execute it"
  echo "# Bash command:  $PYTHON3 $PYPARSER ${pyparser_args[$pi]} -pe ${NPE} -postscript ps${TASK}.sh -psfgen ${CURRPSFGEN} -inpdb ${CURRPDB}"
  $PYTHON3 $PYPARSER ${pyparser_args[$pi]} -pe ${NPE} -postscript $PS -psfgen ${CURRPSFGEN} -inpdb ${CURRPDB}
  #echo "Testing"
  #exit
  
  echo "# Bash command: ./$PS $TASK -task $TASK -nesting-level 2"
  ./$PS -task $TASK -nesting-level 2
  if [ $? -ne 0 ]; then
    echo "Error: Postscript $PS failed; see output immediately above. Exiting."
    exit 1
  fi
  read CURRPSF CURRPDB CURRCFG < .tmpvar
done

# downselect
if [[ "$SEL_STR" != "" ]]; then
  TASK=$((TASK+1))
  echo "# Task $TASK: Downselecting based on selection text $SEL_STR"
  cat > tmp.tcl << EOF
  set selstr [join [split "$SEL_STR" "-"] " "]
  extract_psf_pdb $CURRPSF $CURRPDB \$selstr config${TASK}.psf config${TASK}.pdb
  exit
EOF
  $VMD -dispdev text -e tmp.tcl > ${TASK}-psfgen.log
  CURRPSF=config${TASK}.psf
  CURRPDB=config${TASK}.pdb
  echo tmp.tcl >> .tmpfiles
fi

# solvate
TASK=$((TASK+1))
echo "# Task $TASK: Generating solvated system config${TASK}.psf/.pdb from ${CURRPSF}+${CURRPDB}"
echo "# Bash command: $VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/solv.tcl -args -psf $CURRPSF -pdb $CURRPDB -outpre config${TASK} $CUBICBOX > $solvation_log"
$VMD -dispdev text -e $PSFGEN_BASEDIR/scripts/solv.tcl -args -psf $CURRPSF -pdb $CURRPDB -outpre config${TASK} $CUBICBOX > $solvation_log
if [ $? -ne 0 ]; then
    echo "Error: Solvation using VMD failed.  Check $solvation_log. Exiting."
    exit 1
fi
CURRPSF=config${TASK}.psf
CURRPDB=config${TASK}.pdb

# run solvated NAMD
TASK=$((TASK+1))
if [ ! -f "$CURRCFG" ]; then
    echo "Error: Previous-stage NAMD configuration file $CURRCFG not found."
fi
# get parameter designations from last config file
grep ^parameters $CURRCFG | grep -v water > _par.inp
echo "# No binary inputs yet -- this run begins using PDB coordinates" > _bin.inp
firsttimestep=0
ls=`echo "${#numsteps[@]} - 1" | bc`
for s in `seq 0 $ls`; do
    echo "# Running namd2 (stage $s of ${#numsteps[@]}) on solvated system"
    lastnamd=run${TASK}_stage${s}.namd
    lastsys=config${TASK}_stage${s}
    thislog=run${TASK}_stage${s}.log
    cat $PSFGEN_BASEDIR/templates/solv.namd | \
        sed "/# SYSTEM CONFIGURATION FILES BEGIN/r _bin.inp" | \
        sed "/# SYSTEM CONFIGURATION FILES END/i structure $CURRPSF" | \
        sed "/# SYSTEM CONFIGURATION FILES END/i coordinates $CURRPDB" | \
        sed "/# PARAMETER FILES BEGIN/r _par.inp" | \
        sed s/%STAGE%/${s}/g | \
        sed s/%OUT%/${lastsys}/g | \
        sed s/%NUMSTEPS%/${numsteps[$s]}/g | \
        sed s/%SEED%/${seed}/g | \
        sed s/%TEMPERATURE%/${temperature}/g | \
        sed s/%FIRSTTIMESTEP%/$firsttimestep/g > $lastnamd
    echo "# Bash command:  $CHARMRUN +p${NPE} $NAMD2 $lastnamd > $thislog"
    $CHARMRUN +p${NPE} $NAMD2 $lastnamd > $thislog
    if [ $? -ne 0 ]; then

        echo "Error: NAMD failedat stage $s.  Check log file $thislog. Exiting."
        exit 1
    fi
    firsttimestep=`echo "100 + $firsttimestep + ${numsteps[$s]}" | bc`
    ss=$((s+1))
    
    echo "bincoordinates ${lastsys}.coor" > _bin.inp
    echo "binvelocities  ${lastsys}.vel"  >> _bin.inp
    echo "extendedsystem ${lastsys}.xsc"  >> _bin.inp
done
echo "# All NAMD task$(ess ${#numsteps[@]}) complete.  Packaging for production."

# Prep for production MD
# copy all charmm parameter files to this directory and
# create 'par.inp', which contains local dir file names
${PSFGEN_BASEDIR}/scripts/cp_charmm.sh $lastnamd
firsttimestep=0
cat $PSFGEN_BASEDIR/templates/prod.namd | \
    sed "/# SYSTEM CONFIGURATION FILES BEGIN/r _bin.inp" | \
    sed "/# SYSTEM CONFIGURATION FILES END/i structure $CURRPSF" | \
    sed "/# SYSTEM CONFIGURATION FILES END/i coordinates $CURRPDB" | \
    sed "/# PARAMETER FILES BEGIN/r par.inp" | \
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
echo "# Created namdprod.tgz."
if (( $DO_TOPOGROMACS == 1 )); then
    echo "# Generating Gromacs topology"
    echo "# Bash command: $PSFGEN_BASEDIR/scripts/tg.sh  -psf $CURRPSF -i config${TASK}_stage${s} -top $TG_TOP -pdb $TG_PDB"
    $PSFGEN_BASEDIR/scripts/tg.sh  -psf $CURRPSF -i config${TASK}_stage${s} -top $TG_TOP -pdb $TG_PDB -nesting-level 2
    tar zvcf gmx.tgz $TG_TOP $TG_PDB
    echo "# Created gmx.tgz."
fi
echo "# do_py.sh is done.  Thanks for using do_py.sh."