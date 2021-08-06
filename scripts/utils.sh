function check_command {
    # check for executable named by first argument or contained
    # by a possible environment variable named by the second
    # if the environment variable is not set, set it
    cmd=$1
    envvar=$2
    if ! command -v $cmd &> /dev/null; then
        # no vmd in the users PATH...
        if [[ -z "${!envvar}" ]]; then
            # user did not set a environment variable
            echo "Error: $cmd could not be found"
            exit
        elif [[ ! -f ${!envvar} ]]; then
            echo "No $cmd found at environment variable ${!envvar}"
            exit
        fi
    else
        export ${envvar}=$cmd
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
