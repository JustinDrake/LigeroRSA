#!/bin/bash

# Reading the command line
POSITIONAL=()
FASTTRACK=false
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -f|--fasttrack)
    FASTTRACK=true
    shift                    # past argument
    ;;
    *)                       # unknown option
    POSITIONAL+=("$1")       # save it in an array for later
    shift                    # past argument
esac
done
set -- "${POSITIONAL[@]}"    # restore positional parameters

if [[ $# -ne 2 ]]; then
    echo "USAGE: rebuild_and_run_simulation.sh --fasttrack <key> <simulation_tag>"
    echo
    echo "-f, --fasttrack: builds the binaries and provisions the EC2 infrastructure in parallel"
    exit
fi

if [ "$FASTTRACK" = true ]; then
    echo "Building binaries in the background, see AWS/build.log for status..."
    echo
    export FASTTRACK
    (./rebuild.sh > AWS/build.log 2>&1) &
    ./run_simulation.sh $1 $2

else   
    ./rebuild.sh
    ./run_simulation.sh $1 $2
fi