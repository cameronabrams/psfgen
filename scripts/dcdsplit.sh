#!/bin/bash
# don't use this
if [[ -z "${PSFGEN_BASEDIR}" ]]; then
    PSFGEN_BASEDIR=${HOME}/research/psfgen
    if [[ ! -d $PSFGEN_BASEDIR ]]; then
        echo "Error: No PSFGEN_BASEDIR found."
        exit -1
    fi
fi
source $PSFGEN_BASEDIR/scripts/utils.sh
check_command catdcd
N=1
DCD=""
i=1
ARGC=$#
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-dcd" ]; then
    i=$((i+1))
    DCD=${!i}
  elif [ "${!i}" = "-N" ]; then
    i=$((i+1))
    N=${!i}
  fi
  i=$((i+1))
done
check_file $DCD

Ntot=`catdcd -num $DCD | grep Total | awk '{print $3}'`
Nper=`echo "$Ntot/$N" |bc`
Nrem=`echo "$Ntot-$N*$Nper"|bc`
echo "Ntot $Ntot Nper $Nper Nrem $Nrem"
first=1; # dcd frame numbering begins at "1" not "0"
last=$Nper
slice=1
while [[ $last -le $Ntot ]]; do
   echo "first $first last $last"
   first=$((last+1))
   last=$((last+Nper))
   catdcd -first $first -last $last -o ${DCD%.dcd}-slice${slice}.dcd $DCD
   slice=$((slice+1))
done
if (($Nrem>0)); then
   last=$((first+$Nrem-1))
   echo "rem first $first last $last"
   catdcd -first $first -last $last -o ${DCD%.dcd}-slice${slice}.dcd $DCD
fi

