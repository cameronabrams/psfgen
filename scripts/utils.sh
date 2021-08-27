function check_command {
    # check for executable named by first argument or contained
    # by a possible environment variable named by the second
    # if the environment variable is not set, set it
    envvar="NONE"
    cmd=$1
    if [ "$#" -eq 2 ]; then
        cmd=$1
        envvar=$2
    fi
    if ! command -v $cmd &> /dev/null; then
        if [[ "${envvar}" == "NONE" ]]; then
            # user did not set a environment variable
            echo "Error: $cmd could not be found"
            exit
        else
            #echo "here"
            if [[ ! -f ${!envvar} ]]; then
               echo "No $cmd found at environment variable ${!envvar}"
               exit
            else
               echo "Aliasing $cmd to ${!envvar}"
               alias $cmd=${!envvar}
            fi
        fi
    else
        # if the name of an unset environment variable is provided, set it
        if [[ ! -z "${envvar}" ]]; then
           export ${envvar}=$cmd
        fi
    fi
}

function check_file {
    if [ ! -f $check_file ]; then
       echo "Error: $check_file not found."
       exit 
    fi
}

function ess {
    n=$1
    if (( $n > 1 )); then
        echo "s"
    fi
}

function isare {
    n=$1
    if (( $n > 1 )); then
        echo "are"
    else
        echo "is"
    fi
}

function indent {
    l=$1
    c=$2
    yes "$2" | head -$l | tr -d "\n"
}
