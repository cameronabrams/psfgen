#!/bin/bash

function prep_namd_restart {
    CONF=$1
    LOG=$2
    RECONF=$3
    REOUTNAME=$(basename "$RECONF" | cut -d. -f1)
    if [ ! -f $CONF ] ; then
        echo "error: $CONF not found."
        return 1
    fi
    if [ ! -f $LOG ] ; then
        echo "error: $LOG not found."
        return 1
    fi
    stepsrequested=`grep ^run $CONF | awk '{print $2}' | sed 's/;$//'`
    if [ -z "${stepsrequested}" ]; then
         stepsrequested=`grep ^numsteps $CONF | awk '{print $2}' | sed 's/;$//'`
         if [ -z "${stepsrequested}" ]; then
             echo "Error: $CONF does not contain a run or numsteps statement."
             return 1
         fi
    fi
    lastout=`grep "set outputname" $CONF | awk '{print $3}' | sed 's/;$//'`
    if [ -z "${lastout}" ]; then
         lastout=`grep ^outputname $CONF | awk '{print $2}' | sed 's/;$//'`
         if [ -z "${lastout}" ]; then
             echo "Error: $CONF does not contain an outputname specification."
             return 1
         fi
    fi
    stepsrun=`grep "EXTENDED SYSTEM TO RESTART" $LOG | tail -1 | awk '{print $NF}'`
    if [ -z "${stepsrun}" ]; then
        echo "Error: Cannot determine checkpoint timestep from $LOG"
        return 1
    fi
    stepsleft=$(($stepsrequested-$stepsrun))
    if [[ $stepsleft -eq 0 ]]; then
        echo "${LOG} indicates run has finished ($stepsleft steps left); no restart is necessary."
        return 1
    fi 
    for suf in coor vel xsc; do
        if [ ! -f ${lastout}.restart.${suf} ]; then
            echo "Error: ${lastout}.restart.${suf} not found"
            return 1
        fi
        if [ ! if ${REOUTNAME}.restart.${suf} ]; then
            echo "Error: ${REOUTNAME}.restart.${suf} already exists"
            return 1
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
                sed '4 i bincoordinates '${lastout}'.restart.coor' | \
                sed '5 i binvelocities '${lastout}'.restart.vel' | \
                sed '6 i extendedsystem '${lastout}'.restart.xsc' | \
                sed '/^run/ i firsttimestep '$stepsrun | \
                sed '/^run/ c run '$stepsleft > $RECONF
    echo "Created restart config $RECONF"
    return 0
}


