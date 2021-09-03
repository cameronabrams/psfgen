#!/bin/bash
#
# this program accepts as an argument a NAMD config file, the log from a run generated
# by namd2 run on that config file, and the name of a new config file to create that
# restarts the existing simulation from the last checkpoint, if it did not complete.
#
# some assumptions are made about the syntax of the input config file; read output carefully.
#
# example:
#
# ./prep_namd_restart.sh -i first_run.conf -l run.log -o restart_run.conf [--new-numsteps ###] [--ensemble (nvt/npt)]
#
# Cameron Abrams cfa22@drexel.edu
# 2020-2021

ENSEMBLE="nochange"
i=1
ARGC=$#
quiet=0
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-i" ]; then
    i=$((i+1))
    CONF=${!i}
  elif [ "${!i}" = "-l" ]; then
    i=$((i+1))
    LOG=${!i}
  elif [ "${!i}" = "-o" ]; then
    i=$((i+1))
    RECONF=${!i}
  elif [ "${!i}" = "-q" ]; then
    quiet="1"
  elif [ "${!i}" = "--new-numsteps" ]; then
    i=$((i+1))
    USERSTEPS=${!i}
  elif [ "${!i}" = "--ensemble" ]; then
     i=$((i+1))
     ENSEMBLE=${!i}
  else
    echo "${!i}: not recognized"
  fi
  i=$((i+1))
done

if [[ "$quiet" != "0" ]]; then
    echo "# prep_namd_restart.sh: NAMD continuation utility -- Cameron F Abrams -- cfa22@drexel.edu"
    echo "# Command:"
    echo "#    prep_namd_restart.sh ${@}"
fi

if [[ "$ENSEMBLE" != "nochange" ]] && [[ "$ENSEMBLE" != "nvt" ]] && [[ "$ENSEMBLE" != "npt" ]]; then
    echo "Error: ensemble $ENSEMBLE not recognized"
    exit 1
fi

REOUTNAME=$(basename "$RECONF" | cut -d. -f1)
if [ ! -f $CONF ] ; then
    echo "Error: $CONF not found."
    exit 1
fi
if [ ! -f $LOG ] ; then
    echo "Error: $LOG not found."
    exit 1
fi
HAS_LANGEVIN_THERMOSTAT=`grep -c "LANGEVIN DYNAMICS ACTIVE" $LOG`
HAS_LANGEVIN_BAROSTAT=`grep -c "LANGEVIN PISTON PRESSURE CONTROL ACTIVE" $LOG`
CURRENT_ENSEMBLE="unknown"
if [[ "$HAS_LANGEVIN_THERMOSTAT" == "1" ]]; then
    CURRENT_ENSEMBLE="nvt"
    if [[ "$HAS_LANGEVIN_BAROSTAT" == "1" ]]; then
        CURRENT_ENSEMBLE="npt"
    fi
fi
if [[ "$ENSEMBLE" != "nochange" ]] && [[ "$ENSEMBLE" != "$CURRENT_ENSEMBLE" ]]; then
   echo "Ensemble change from $CURRENT_ENSEMBLE to $ENSEMBLE requested"
fi
#echo "CURRENT_ENSEMBLE: $CURRENT_ENSEMBLE"
#exit
stepsrun=`grep "EXTENDED SYSTEM TO RESTART" $LOG | tail -1 | awk '{print $NF}'`
if [ -z "${stepsrun}" ]; then
    echo "Error: Cannot determine checkpoint timestep from $LOG"
    exit 1
fi
if [ ! -z "${USERSTEPS}" ]; then
    stepsleft=$USERSTEPS
else
    stepsrequested=`grep ^run $CONF | awk '{print $2}' | sed 's/;$//'`
    if [ -z "${stepsrequested}" ]; then
        stepsrequested=`grep ^numsteps $CONF | awk '{print $2}' | sed 's/;$//'`
        if [ -z "${stepsrequested}" ]; then
            echo "Error: $CONF does not contain a run or numsteps statement."
            exit 1
        fi
    else
        # if there is a run statement, the total number of steps would be the
        # number specified in the run statement PLUS the value of firsttimestep
        firsttimestep=`grep ^firsttimestep $CONF | awk '{print $2}' | sed 's/;$//'`
        if [ -z "${firsttimestep}" ]; then
            firsttimestep=0
        else
            echo "$CONF contains a firsttimestep $firsttimestep and a run $stepsrequested"
        fi
        echo "firsttimestep $firsttimestep"
        stepsrequested=$(($stepsrequested+$firsttimestep))
    fi
    stepsleft=$(($stepsrequested-$stepsrun))
fi
lastout=`grep "set outputname" $CONF | awk '{print $3}' | sed 's/;$//'`
if [ -z "${lastout}" ]; then
    lastout=`grep ^outputname $CONF | awk '{print $2}' | sed 's/;$//'`
    if [ -z "${lastout}" ]; then
        echo "Error: $CONF does not contain an outputname specification."
        exit 1
    fi
fi
if [[ $stepsleft -eq 0 ]]; then
    echo "${LOG} indicates run has finished ($stepsleft steps left); no restart is necessary."
    echo "The final checkpoint is in ${lastout}.coor, ${lastout}.vel, and ${lastout}.xsc."
    exit 1
fi 
tmdon=`grep "^tmd on" $CONF | awk '{print $2}'`
if [ ! -z "${tmdon}" ]; then
    # check config for a tmdinitialrmsd
    tmdinitialrmsd=`grep ^tmdinitialrmsd $CONF | awk '{print $2}'`
    if [ -z "${tmdinitialrmsd}" ]; then
        # try to find it in the log
        tmdinitialrmsd=`grep ^TMD $LOG | grep "Domain: 0" | head -1 | awk '{print $5}'`
        if [ -z "${tmdinitialrmsd}" ]; then
            echo "Error: $CONF indicates a TMD run but neither $CONF nor $LOG indicates tmdinitialrmsd."
            exit 1
        fi
    else
        tmdinitialrmsd_inconf=1
    fi
fi

RESINFILENAME=""
for suf in coor vel xsc; do
    if [ -f ${lastout}.restart.${suf} ]; then
        RESINFILENAME=".restart"
    elif [ -f ${lastout}.${suf} ]; then
        RESINFILENAME=""
    else
        echo "Error: Last checkpoint file ${lastout}.restart.${suf} or ${lastout}.${suf} not found."
        exit 1
    fi
    if [ -f ${REOUTNAME}.restart.${suf} ]; then
        echo "Error: Next checkpoint ${REOUTNAME}.restart.${suf} already exists."
        exit 1
    fi
done
cat $CONF | sed '/^#/d' | sed '/^$/d' | sed '/^firsttimestep/d' | \
            sed '/^temperature/d' | \
            sed '/source cell/d' | \
            sed '/^bincoordinates/d' | \
            sed '/^binvelocities/d' | \
            sed '/^extendedsystem/d' | \
            sed '1 i # restart file generated from '$CONF' and '$LOG | \
            sed '/set outputname/ c set outputname '$REOUTNAME | \
            sed '4 i bincoordinates '${lastout}${RESINFILENAME}'.coor' | \
            sed '5 i binvelocities '${lastout}${RESINFILENAME}'.vel' | \
            sed '6 i extendedsystem '${lastout}${RESINFILENAME}'.xsc' | \
            sed '/^run/ i firsttimestep '$stepsrun | \
            sed '/^run/ c run '$stepsleft > $RECONF
if [[ $CURRENT_ENSEMBLE == "npt" ]] && [[ $ENSEMBLE == "nvt" ]]; then
    cat $RECONF | sed 's/^\<langevinpiston\>.*/langevinpiston off/' > tmp
    mv tmp $RECONF
fi
if [ ! -z "${tmdon}" ]; then
    if [ -z "${tmdinitialrmsd_inconf}" ]; then
        cat $RECONF | sed '/tmd on/ a tmdinitialrmsd '${tmdinitialrmsd} > tmp
        mv tmp $RECONF
    fi
fi	
echo "Created restart config $RECONF"
