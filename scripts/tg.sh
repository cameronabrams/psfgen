#!/bin/bash
# topogromacs control script
# cameron f abrams cfa22@drexel.edu
if [[ -z "${PSFGEN_BASEDIR}" ]]; then
    PSFGEN_BASEDIR=${HOME}/research/psfgen
    if [[ ! -d $PSFGEN_BASEDIR ]]; then
        echo "Error: No PSFGEN_BASEDIR found."
        exit -1
    fi
fi
source $PSFGEN_BASEDIR/scripts/utils.sh
PDB=""
TOP="none"
INPUTNAME="none"
PSF="none"
MDP=""
TPR=""
LOG="topogromacs.log"
INTERPDB="tg_needsbox.pdb"
CELLDIMFILE="tg-cell-nm.in"
SYSTEMNAME="cfa-tg-psfgen"
nesting_level=1
i=1
ARGC=$#
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-pdb" ]; then
    i=$((i+1))
    PDB=${!i}
  elif [ "${!i}" = "-psf" ]; then
    i=$((i+1))
    PSF=${!i}
  elif [ "${!i}" = "-top" ]; then
    i=$((i+1))
    TOP=${!i}
  elif [ "${!i}" = "-mdp" ]; then
    i=$((i+1))
    MDP=${!i}
  elif [ "${!i}" = "-tpr" ]; then
    i=$((i+1))
    TPR=${!i}
  elif [ "${!i}" = "-log" ]; then
    i=$((i+1))
    LOG=${!i}
  elif [ "${!i}" = "--cell-dim-file" ]; then
    i=$((i+1))
    CELLDIMFILE=${!i}
  elif [ "${!i}" = "--interpdb" ]; then
    i=$((i+1))
    INTERPDB=${!i}
  elif [ "${!i}" = "--systemname" ]; then
    i=$((i+1))
    SYSTEMNAME=${!i}
  elif [ "${!i}" = "-i" ]; then
    i=$((i+1))
    INPUTNAME=${!i}
  elif [ "${!i}" = "-nesting-level" ]; then
    i=$((i+1))
    nesting_level=${!i}
  else
    echo "${!i} unknown."
  fi
  i=$((i+1))
done

check_command vmd
check_command gmx

INDENT=`indent $nesting_level "#"`

if [ "$PSF" = "none" ]; then
   echo "Error: You must specify the input psf file name with the -psf option."
   exit
fi
if [ "$INPUTNAME" = "none" ]; then
   echo "Error: You must specify the input name (prefix of coor and xsc files) with the -i option."
   exit
fi
if [ "$TOP" = "none" ]; then
   echo "Error: You must specify the output gromacs topology file name with the -top option."
   exit
fi
if [ "$PDB" = "none" ]; then
   echo "Error: You must specify the output gromacs pdb file name with the -pdb option."
   exit
fi
for f in $PSF ${INPUTNAME}.coor ${INPUTNAME}.xsc; do
  if [ ! -f $f ]; then
    echo "Error: $f not found.  Cannot execute topogromacs."
    exit
  fi
done

echo "Executing topogromacs VMD script to convert $PSF/$INPUTNAME to $TOP/$INTERPDB"
echo "Command: vmd -dispdev text -e $PSFGEN_BASEDIR/scripts/tg.tcl -args -psf $PSF -top $TOP -i $INPUTNAME -opdb $INTERPDB --cell-dim-file $CELLDIMFILE"
vmd -dispdev text -e $PSFGEN_BASEDIR/scripts/tg.tcl -args -psf $PSF -top $TOP -i $INPUTNAME -opdb $INTERPDB --cell-dim-file $CELLDIMFILE
if [ $? -ne 0 ]; then
    echo "Topogromacs script failed."
    exit 1
fi
echo "Calling gmx editconf to combine box size info from $CELLDIMFILE with $INTERPDB to generate $PDB"
echo "Command: gmx editconf -f $INTERPDB -o $PDB -box `cat $CELLDIMFILE` 2>&1 >> $LOG"
gmx editconf -quiet -f $INTERPDB -o $PDB -box `cat $CELLDIMFILE` 2>&1 >> $LOG
echo "Done.  Results in $LOG."
if [ "$MDP" != "" ] && [ "$TPR" != "" ] && [ $TOP != "" ]; then
  for f in $MDP $PDB $TOP; do
    if [ ! -f $f ]; then
      echo "Error: $f not found.  Cannot execute gmx grompp."
      exit
    fi
  done
  cat $TOP | sed s/vmdmolecule2/$SYSTEMNAME/ > .tmp; mv .tmp $TOP
  gmx grompp -quiet -f $MDP -c $PDB -p $TOP -o $TPR -maxwarn 2
else
   echo "Next command: gmx grompp -f whatever.mdp -c $PDB -p $TOP -o whatever.tpr -maxwarn 2"
fi

