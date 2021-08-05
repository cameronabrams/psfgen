#!/bin/bash
# topogromacs control script
# cameron f abrams cfa22@drexel.edu

PDB="tg_coor.pdb"
TOP="tg_top.pdb"
INPUTNAME="none"
PSF="none"
MDP=""
TPR=""
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

echo "Executing topogromacs VMD script to convert $PSF/$INPUTNAME to $TOP/$PDB..."
vmd -dispdev text -e $PSFGEN_BASEDIR/scripts/tg.tcl -args -psf $PSF -top $TOP -i $INPUTNAME 2>&1 > topogromacs.log
if [ $? -ne 0 ]; then
    echo "Topogromacs script failed. Check topogromacs.log."
    exit 1
fi
gmx editconf -f tg_needsbox.pdb -o $PDB -box `cat tg-cell-nm.in` 2>&1 >> topogromacs.log
echo "Done.  Results in topogromacs.log."
echo "Do this next: gmx grompp -f whatever.mdp -c $PDB -p $TOP -o whatever.tpr -maxwarn 2"
# example grompp for production md
#gmx grompp -f $MDP -c $PDB -p $TOP -o $TPR -maxwarn 2
