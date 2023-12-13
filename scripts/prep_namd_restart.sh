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

banner () {
    echo "# prep_namd_restart.sh: NAMD continuation utility -- Cameron F Abrams -- cfa22@drexel.edu"
    echo "# Command:"
    echo "#    prep_namd_restart.sh ${@}"
}

help () {
    echo "# NAMD continuation utility -- Cameron F Abrams -- cfa22@drexel.edu"
    echo "Usage:"
    echo " $ prep_namd_restart.sh -i <existing_namd_config> -l <existing_namd_log> -o <new_namd_config_to_generate> [options]"
    echo "Options:"
    echo "    --new-numsteps <#>             change the number of timesteps to this value"
    echo "    --addsteps <#>                 change the number of timesteps by adding this value to the current number"
    echo "    --ensemble <nochange|nvt|npt>  set the ensemble; default is nochange"
    echo "    -q                             supress messages"
    echo "    --help                         generate this message"
}

ENSEMBLE="nochange"
USERADDSTEPS=0
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
  elif [ "${!i}" = "--addsteps" ]; then
    i=$((i+1))
    USERADDSTEPS=${!i}
  elif [ "${!i}" = "--ensemble" ]; then
    i=$((i+1))
    ENSEMBLE=${!i}
  elif [ "${!i}" = "--help" ]; then
    help
    exit
  else
    echo "${!i}: not recognized"
  fi
  i=$((i+1))
done

cfgext=`echo $RECONF| awk -F. '{print $NF}'`
REOUTNAME=`echo $RECONF | sed s/".${cfgext}//"`
FILES=($RECONF)

for pfile in `grep -i ^parameters $CONF | awk '{print $NF}'`; do
   FILES+=($pfile)
done
FILES+=(`grep -i ^structure $CONF|awk '{print $NF}'`)
FILES+=(`grep -i ^coordinates $CONF|awk '{print $NF}'`)
FILES+=(`grep -i ^colvarsconfig $CONF|awk '{print $NF}'`)
FILES+=(`grep -iw ^tmdfile $CONF|awk '{print $NF}'`)
FILES+=(`grep -iw ^tmdfile2 $CONF|awk '{print $NF}'`)

if [[ "$quiet" != "0" ]]; then
    banner
fi

if [[ "$ENSEMBLE" != "nochange" ]] && [[ "$ENSEMBLE" != "nvt" ]] && [[ "$ENSEMBLE" != "npt" ]]; then
    echo "Error: ensemble $ENSEMBLE not recognized"
    exit 1
fi

if [ ! -f $CONF ] ; then
    echo "Error: $CONF not found."
    exit 1
fi
if [ ! -f $LOG ] ; then
    echo "Error: $LOG not found."
    exit 1
fi
HAS_LANGEVIN_THERMOSTAT=`grep -c "LANGEVIN DYNAMICS ACTIVE" $LOG`
# Info: LANGEVIN TEMPERATURE   310
LANGEVIN_TEMPERATURE=`grep "Info: LANGEVIN TEMPERATURE" $LOG | awk '{print $NF}'`
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
            echo "Previous $CONF contains a firsttimestep $firsttimestep and a run $stepsrequested"
        fi
        #echo "firsttimestep $firsttimestep"
        stepsrequested=$(($stepsrequested+$firsttimestep))
    fi
    stepsleft=$(($stepsrequested-$stepsrun))
    stepsleft=$(($stepsleft+$USERADDSTEPS))
fi
lastout=`grep "set outputname" $CONF | grep -v \# | awk '{print $3}' | sed 's/;$//'`
if [ -z "${lastout}" ]; then
    lastout=`grep -i ^outputname $CONF | awk '{print $2}' | sed 's/;$//'`
    if [ -z "${lastout}" ]; then
        echo "Error: $CONF does not contain an outputname specification."
        exit 1
    fi
fi
if [[ $stepsleft -eq 0 ]]; then
    echo "${LOG} indicates run has finished ($stepsleft out of $stepsrequested steps left); no restart is necessary."
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
REINITVELS="no"
for suf in coor vel xsc; do
    if [ -f ${lastout}.restart.${suf} ]; then
        RESINFILENAME=".restart"
    elif [ -f ${lastout}.${suf} ]; then
        RESINFILENAME=""
    else
        echo "Last checkpoint file ${lastout}.restart.${suf} or ${lastout}.${suf} not found."
        if [[ $suf == "vel" ]]; then
           REINITVELS="yes"
        else
           echo "Error."
           exit 1
        fi
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
            sed '/^extendedSystem/d' | \
            sed '1 i # restart file generated from                    associates={}
 '$CONF' and '$LOG | \
            sed '/set outputname/ c set outputname '$REOUTNAME | \
            sed '/^outputName/ c outputName \$outputname' | \
            sed '/^outputname/ c outputname \$outputname' | \
            sed '4 i bincoordinates '${lastout}${RESINFILENAME}'.coor' | \
            sed '5 i extendedsystem '${lastout}${RESINFILENAME}'.xsc' | \
            sed '6 i set outputname '$REOUTNAME | \
            sed '/^run/ i firsttimestep '$stepsrun | \
            sed '/^run/ c run '$stepsleft | \
            sed '/^numsteps/ c run '$stepsleft > $RECONF
if [[ $CURRENT_ENSEMBLE == "npt" ]] && [[ $ENSEMBLE == "nvt" ]]; then
    cat $RECONF | sed 's/^\<langevinpiston\>.*/langevinpiston off/' > tmp
    mv tmp $RECONF
fi
if [[ $REINITVELS == "yes" ]]; then 
    cat $RECONF | sed '3 i temperature '$LANGEVIN_TEMPERATURE > tmp
    mv tmp $RECONF
else
    cat $RECONF | sed '5 i binvelocities '${lastout}${RESINFILENAME}'.vel' > tmp
    mv tmp $RECONF
fi
if [ ! -z "${tmdon}" ]; then
    if [ -z "${tmdinitialrmsd_inconf}" ]; then
        cat $RECONF | sed '/tmd on/ a tmdinitialrmsd '${tmdinitialrmsd} > tmp
        mv tmp $RECONF
    fi
fi
FILES+=(`grep -i ^bincoordinates $RECONF|awk '{print $NF}'`)
FILES+=(`grep -i ^binvelocities $RECONF|awk '{print $NF}'`)
FILES+=(`grep -i ^extendedsystem $RECONF|awk '{print $NF}'`)

echo ${FILES[@]} > ${RECONF}.files
echo "Created restart config $RECONF"
