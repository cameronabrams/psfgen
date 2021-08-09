#!/bin/bash
# topogromacs control script
# cameron f abrams cfa22@drexel.edu

PDB=""
TOP="none"
INPUTNAME="none"
PSF="none"
MDP=""
TPR=""
LOG="topogromacs.log"
i=1
ARGC=$#
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-pdb" ]; then
    i=$((i+1))
    PDB=${!i}
  fi
  if [ "${!i}" = "-psf" ]; then
    i=$((i+1))
    PSF=${!i}
  fi
  if [ "${!i}" = "-top" ]; then
    i=$((i+1))
    TOP=${!i}
  fi
  if [ "${!i}" = "-mdp" ]; then
    i=$((i+1))
    MDP=${!i}
  fi
  if [ "${!i}" = "-tpr" ]; then
    i=$((i+1))
    TPR=${!i}
  fi
  if [ "${!i}" = "-log" ]; then
    i=$((i+1))
    LOG=${!i}
  fi
  if [ "${!i}" = "-i" ]; then
    i=$((i+1))
    INPUTNAME=${!i}
  fi
  i=$((i+1))
done

for cmd in vmd gmx ; do
    if ! command -v $cmd &> /dev/null
    then
        echo "Error: $cmd could not be found"
        exit
    fi
done

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

echo "Executing topogromacs VMD script to convert $PSF/$INPUTNAME to $TOP/$PDB..."
vmd -dispdev text -e $PSFGEN_BASEDIR/scripts/tg.tcl -args -psf $PSF -top $TOP -i $INPUTNAME 2>&1 > $LOG
if [ $? -ne 0 ]; then
    echo "Topogromacs script failed. Check topogromacs.log."
    exit 1
fi
gmx editconf -f tg_needsbox.pdb -o $PDB -box `cat tg-cell-nm.in` 2>&1 >> $LOG
echo "Done.  Results in $LOG."
if [ "$MDP" != "" ] && [ "$TPR" != "" ] && [ $TOP != "" ]; then
  for f in $MDP $PDB $TOP; do
    if [ ! -f $f ]; then
      echo "Error: $f not found.  Cannot execute gmx grompp."
      exit
    fi
  done
  gmx grompp -f $MDP -c $PDB -p $TOP -o $TPR -maxwarn 2
else
   echo "Do this next: gmx grompp -f whatever.mdp -c $PDB -p $TOP -o whatever.tpr -maxwarn 2"
fi

