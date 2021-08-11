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
if [[ -z "${PSFGEN_BASEDIR}" ]]; then
    PSFGEN_BASEDIR=${HOME}/research/psfgen
    if [[ ! -d $PSFGEN_BASEDIR ]]; then
        echo "Error: No PSFGEN_BASEDIR found."
        exit -1
    fi
fi
source $PSFGEN_BASEDIR/scripts/utils.sh

PDB="tg_coor.pdb"
TPR="tg.tpr"
TRR="tg.trr"
PSF="tg_redordered.psf"
CSELSTR=""
DCD=""
APPEND=1

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
  if [ "${!i}" = "-no-append" ]; then
    APPEND=0
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

check_command vmd
check_command gmx
check_command catdcd

gmx dump -s $TPR -om tmp.mdp
dt=`grep ^dt tmp.mdp|awk '{print $3}'`
nstxout=`grep -w ^nstxout tmp.mdp| grep -v compressed|awk '{print $3}'`
nsteps=`grep ^nsteps tmp.mdp|awk '{print $3}'`
#echo "$dt $nsteps $nstxout ($nsteps/$nstxout)+1"
nframes_expected=`echo "($nsteps/$nstxout)+1"|bc`
frame_interval_ps=`echo "$dt*$nstxout"|bc`
echo "$TPR: $nframes_expected frames are in $TRR, interval $frame_interval_ps ps"
b=0

if [ -f $DCD ]; then
   if $APPEND; then
    nframes_converted=`catdcd -num $DCD | grep Total | awk '{print $3}'`
    echo "# Output $DCD exists with $nframes_converted frames."
    b=`echo "$nframes_converted*$nstxout*$dt" | bc -l`
    echo "# Will begin reading $TRR at time $b ps."
    echo "# Copying $DCD to PREV-$DCD..."
    cp $DCD PREV-${DCD}
  else
    echo "$DCD found; will be overwritten"
  fi
fi

# make appropriate ndx files
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

# center on the centersel
gmx trjconv -f $TRR -b $b -s $PDB -o tmp1.trr -pbc nojump -center -n all.ndx < tmp
cat > tmp << EOF
0
EOF
# put molecules back together
gmx trjconv -f tmp1.trr -b $b -s $TPR -o tmp2.trr -pbc mol < tmp

# convert to dcd
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
if [ -f PREV-${DCD} ] && $APPEND; then
    catdcd -o tmp.dcd PREV-${DCD} ${DCD}
    mv tmp.dcd ${DCD}
fi
rm -f tmp tmp.tcl tmp1.trr tmp2.trr

