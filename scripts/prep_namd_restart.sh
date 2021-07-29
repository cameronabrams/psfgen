#!/bin/bash
#
# this function accepts as an argument a NAMD config file, the log from a run generated
# by namd2 run on that config file, and the name of a new config file to create that
# restarts the existing simulation from the last checkpoint, if it did not complete.
#
# some assumptions are made about the syntax of the input config file; read output carefully.
#
# ***To use this function, source this file in your current bash session, or in your .bashrc***
#
# example:
#
# prep_namd_restart first_run.conf run.log restart_run.conf
#
# Cameron Abrams cfa22@drexel.edu
# 2020

function prep_namd_restart {
    CONF=$1
    LOG=$2
    RECONF=$3
    USERSTEPS=$4
    REOUTNAME=$(basename "$RECONF" | cut -d. -f1)
    if [ ! -f $CONF ] ; then
        echo "Error: $CONF not found."
        return 1
    fi
    if [ ! -f $LOG ] ; then
        echo "Error: $LOG not found."
        return 1
    fi
    stepsrun=`grep "EXTENDED SYSTEM TO RESTART" $LOG | tail -1 | awk '{print $NF}'`
    if [ -z "${stepsrun}" ]; then
        echo "Error: Cannot determine checkpoint timestep from $LOG"
        return 1
    fi
    if [ ! -z "${USERSTEPS}" ]; then
        stepsleft=$USERSTEPS
    else
        stepsrequested=`grep ^run $CONF | awk '{print $2}' | sed 's/;$//'`
        if [ -z "${stepsrequested}" ]; then
            stepsrequested=`grep ^numsteps $CONF | awk '{print $2}' | sed 's/;$//'`
            if [ -z "${stepsrequested}" ]; then
                echo "Error: $CONF does not contain a run or numsteps statement."
                return 1
            fi
        else
            # if there is a run statement, the total number of steps would be the
            # number specified in the run statement PLUS the value of firsttimestep
            firsttimestep=`grep ^firsttimestep $CONF | awk '{print $2}' | sed 's/;$//'`
            if [ -z "${firstttimestep}" ]; then
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
            return 1
        fi
    fi
    if [[ $stepsleft -eq 0 ]]; then
        echo "${LOG} indicates run has finished ($stepsleft steps left); no restart is necessary."
        echo "The final checkpoint is in ${lastout}.coor, ${lastout}.vel, and ${lastout}.xsc."
        return 1
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
	            return 1
	        fi
        else
	        tmdinitialrmsd_inconf=1
	    fi
    fi

    for suf in coor vel xsc; do
        if [ ! -f ${lastout}.restart.${suf} ]; then
            echo "Error: Last checkpoint file ${lastout}.restart.${suf} not found."
            return 1
        fi
        if [ -f ${REOUTNAME}.restart.${suf} ]; then
            echo "Error: Next checkpoint ${REOUTNAME}.restart.${suf} already exists."
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
    if [ ! -z "${tmdon}" ]; then
        if [ -z "${tmdinitialrmsd_inconf}" ]; then
            cat $RECONF | sed '/tmd on/ a tmdinitialrmsd '${tmdinitialrmsd} > tmp
	        mv tmp $RECONF
	    fi
    fi	
    echo "Created restart config $RECONF"
    return 0
}


