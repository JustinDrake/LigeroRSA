#!/bin/bash

if [ $# -lt 2 ]; then
    echo "USAGE: run_EC2 <instance> <command> <opt:force foreground>"
    exit 1
fi

INSTANCE=$(cat AWS/Registered_Instances/$1)

if [ $# -eq 3 ]; then
 ssh -o StrictHostKeyChecking=no -i "AWS/$EC2_KEY_ID.pem" $INSTANCE "$2"
else
 ssh -o StrictHostKeyChecking=no -f -i "AWS/$EC2_KEY_ID.pem" $INSTANCE "$2"
fi
