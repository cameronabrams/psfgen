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
for f in `grep ^parameters $CONF | awk -F\) '{print $2}'`; do
  cp ${HOME}${f} .
  ff=`echo $f | awk -F'/' '{print $NF}'`
  echo "parameters $ff" >> par.inp
done

