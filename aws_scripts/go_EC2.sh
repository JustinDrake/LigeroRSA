#!/bin/bash

if [ $# -ne 1 ]; then
    echo "USAGE: go_EC2 <reg_nb>"
    exit 1
fi

INSTANCE=$(cat AWS/Registered_Instances/$1)
echo $INSTANCE
ssh -o StrictHostKeyChecking=no -i "AWS/$EC2_KEY_ID.pem" $INSTANCE
