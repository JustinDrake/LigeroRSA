#!/bin/bash

# Test arguments are properly formed
if [ $# -ne 2 ]; then
    echo "USAGE: run_simulation.sh <key:string> <tag:string>"
    exit 1
fi

# Acquire all parameters
source AWS/config $1 $2

if [ ! "$FASTTRACK" = true ]; then
    if [ ! -f "$BUILD_DIR/$COORDINATOR_BINARY" ] || [ ! -f "$BUILD_DIR/$PARTY_BINARY" ]; then
            echo "Error: binaries missing, aborting now."
            exit
    fi
fi

if [ -f "AWS/$EC2_KEY_ID.pem" ]; then

    # Launch the sim
    echo
    echo "Provisioning the AWS infrastructure"
    echo "-----------------------------------"
    ./rsa_ceremony_deploy.sh
    
    echo
    echo "Starting the simulation"
    echo "-----------------------"
    ./rsa_ceremony_start.sh
    
    echo
    echo "Simulation completed, cleaning up (see AWS/cleanup.log)"
    echo "-------------------------------------------------------"
    ./rsa_terminate_all.sh > AWS/cleanup.log
    
    echo
    echo "All set, check S3 bucket ligero-simulations for analytics"
    echo "---------------------------------------------------------"

else 
    echo "Error: AWS private key missing, aborting now."
fi