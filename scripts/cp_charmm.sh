#!/bin/bash
# 
# cp any system parameter files in the namd script to the current directory
#
# cameron f abrams cfa22@drexel.edu
#
CONF=$1
if [ -f par.inp ]; then
   rm par.inp
fi
touch par.inp
for f in `grep ^parameters $CONF | awk '{print $2}'`; do
  ev=`echo $f | awk -F'/' '{print $1}'`
  if [[ $ev == '$env(HOME)' ]]; then
      pre=${HOME}
  fi
  if [[ $ev == '$env(PSFGEN_BASEDIR)' ]]; then
      pre=${PSFGEN_BASEDIR}
  fi  
  suf=`echo $f | cut -d/ -f 2-`
  cp $pre/$suf .
  ff=`echo $suf | awk -F'/' '{print $NF}'`
  echo "parameters $ff" >> par.inp
done

