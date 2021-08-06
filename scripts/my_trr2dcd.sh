#!/bin/bash

# A tool for converting trr trajectory files from gromacs to dcd files for vmd
# Allows for centering on a custom user selection; selection string must be gromacs-format
#
# Example: Say you have abc.trr generated from abc.pdb and abc.tpr, and you have a congruent psf file abc.ps
# (which is possible if this system was generated from a NAMD system using topogromacs)
# Then you can generate abc.dcd via
#
# > $PSFGEN_BASEDIR/scripts/my_trr2dcd.sh -pdb abc.pdb -trr abc.trr -tpr abc.tpr -psf abc.psf -dcd abc.dcd
#
# Cameron F Abrams cfa22@drexel.edu
# 2021

PDB="tg_coor.pdb"
TPR="tg.tpr"
TRR="tg.trr"
PSF="tg_redordered.psf"
CSELSTR=""
DCD=""

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
  if [ "${!i}" = "-trr" ]; then
    i=$((i+1))
    TRR=${!i}
  fi
  if [ "${!i}" = "-tpr" ]; then
    i=$((i+1))
    TPR=${!i}
  fi
  if [ "${!i}" = "-center-sel-str" ]; then
    i=$((i+1))
    CSELSTR=${!i}
  fi
  if [ "${!i}" = "-dcd" ]; then
    i=$((i+1))
    DCD=${!i}
  fi
  i=$((i+1))
done

if [ "$DCD" = "" ]; then
   DCD=${TRR%trr}.dcd
fi

for infile in $PSF $TRR $PDB $TPR; do
   if [ ! -f $infile ]; then
       echo "Error: $infile not found."
       exit
    fi
done

for cmd in vmd gmx ; do
    if ! command -v $cmd &> /dev/null
    then
        echo "Error: $cmd command could not be found."
        exit
    fi
done
echo "q" > tmp
gmx make_ndx -f $PDB -o tst.ndx < tmp
if [ "$CSELSTR" != "" ]; then
   gmx select -s $PDB -select "$CSELSTR" -on o.ndx
   cat o.ndx tst.ndx > all.ndx
else
   mv tst.ndx all.ndx
fi
cat > tmp << EOF
1
0
EOF
gmx trjconv -f $TRR -s $PDB -o tmp1.trr -pbc nojump -center -n all.ndx < tmp
cat > tmp << EOF
0
EOF
gmx trjconv -f tmp1.trr -s $TPR -o tmp2.trr -pbc mol < tmp
cat > tmp.tcl << EOF
mol new $PSF
mol addfile tmp2.trr waitfor all
set n0 [molinfo top get numframes]
set a [atomselect top all]
set n [animate write dcd $DCD beg 0 end [expr \$n0 - 1] sel \$a top]
puts "\$n frames written to $DCD"
exit
EOF
vmd -dispdev text -e tmp.tcl
rm -f tmp tmp.tcl tmp1.trr tmp2.trr

