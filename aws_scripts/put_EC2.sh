#!/bin/bash

if [ $# -ne 3 ]; then
    echo "USAGE: put_EC2 <reg_nb> <path_here> <path_server>"
    exit 1
fi

INSTANCE=$(cat AWS/Registered_Instances/$1)
rsync -e "ssh -o StrictHostKeyChecking=no -i AWS/$EC2_KEY_ID.pem -q" $2 $INSTANCE:/home/ubuntu/$3 > /dev/null
