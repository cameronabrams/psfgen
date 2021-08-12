#!/bin/bash
#
# MMGBSA_DOIT  (c) 2019 cameron f abrams cfa22@drexel.edu
#
# Use namd2 to peform MMBGSA interaction energy calculations
# on DCD trajectories extracted from raw simulation output DCD's
# using STRIPSPLIT.  This operates on a single replica in a single
# system.

if [[ -z "${PSFGEN_BASEDIR}" ]]; then
    PSFGEN_BASEDIR=${HOME}/research/psfgen
    if [[ ! -d $PSFGEN_BASEDIR ]]; then
        echo "Error: No PSFGEN_BASEDIR found."
        exit -1
    fi
fi
source $PSFGEN_BASEDIR/scripts/utils.sh

# these two commands must be in your PATH
check_command vmd
check_command namd2

TEMPLATECF="$PSFGEN_BASEDIR/scripts/mmgbsa/mmgbsa_template.namd"
DOCALC="YES"
FORCE="YES"
final_results_file="results.rae"
Asel=""
Bsel=""
Aname="A"
Bname="B"
ABname="AB"
NPE=1
stride=1

while [ "$#" -gt 0 ]; do
  case "$1" in
    -psf) PSF="$2"; shift 2;;
    -dcd) DCD="$2"; shift 2;;
    -Asel) Asel="$2"; shift 2;;
    -Bsel) Bsel="$2"; shift 2;;
    -Aname) Aname="$2"; shift 2;;
    -Bname) Bname="$2"; shift 2;;
    -ABname) ABname="$2"; shift 2;;
    -NPE) NPE="$2"; shift 2;;
    -stride) stride="$2"; shift 2;;
    -o) final_results_file="$2"; shift 2;;
    --force) FORCE="YES"; shift 1;;
    -namd-config-template) TEMPLATECF="$2"; shift 2;;
    *) echo "unrecognized argument: $1"
  esac
done

# check to see if calculation was already performed for this replica
if [ -f $final_results_file ]; then
    cp $final_results_file ${final_results_file}.bak
    echo "$final_results_file copied to ${final_results_file}.bak"
    if [ "$FORCE" == "NO" ]; then 
      echo "Final results $final_results_file already exists.  Use --force to force a recalculation."
      DOCALC=NO
    else
      echo "Recalculating."
    fi
fi

# stripsplit
vmd -dispdev text -e $PSFGEN_BASEDIR/scripts/mmgbsa/stripsplit.tcl -args -Asel $Asel -Bsel $Bsel -psf $PSF -stride $stride -Aname $Aname -Bname $Bname -ABname $ABname $DCD

# generate the config file for each type of system, run namd2 to compute energies on existing trajectory, 
# extract potential energy from ENERGY lines in namd2 log
if [ "$DOCALC" == "YES" ]; then 
  for sys in $Aname $Bname $ABname; do
    pf=m-${sys}
    c=${pf}.namd
    l=${pf}.log
    e=${pf}.e
    cat ${TEMPLATECF} | sed s/%SYS%/${sys}/g  > $c
    echo "Running $namd2 +p${NPE} $c"
    namd2 +p${NPE} $c > $l
    echo "Generated $l"
    if [ "$sys" == "$Aname" ] ; then
      grep ^ENERGY $l | awk '{print $2,$14}' > $e
    else
      grep ^ENERGY $l | awk '{print $14}' > $e
    fi
  done
fi

# perform the running average of the difference (complex)-((ligand)+(target)) 
# potential energies
paste m-${Aname}.e m-${Bname}.e m-${ABname}.e | \
      awk 'BEGIN{ra=0.0} {ra+=($4-$3-$2); print $1,ra/NR}' \
      > $final_results_file
echo "Generated $final_results_file."

