#!/bin/bash

if [ $# -lt 1 ]; then

    NB=0
    while IFS=, read -r location location_description instance_type nb_instances ami security_group
        do    
                NB=$((NB+$nb_instances))
        done < AWS/map_deployment_rsa_ceremony
    NB=$((NB-1))
    echo "No number of parties specified, will use $NB parties based on configuration file"
 
else 
    NB=$1
fi

PRIVATE_IP=$(cat AWS/Registered_Private_IPs/1)

echo $COORDINATOR_BINARY
echo $PRIVATE_IP
echo $NB

# Execution
PID=$(./run_EC2.sh 1 "source config \"NULL\" $TAG && ulimit -n 4096 && ulimit -s unlimited && tmux new-session -d -s "rsamain" ./$COORDINATOR_BINARY --ip $PRIVATE_IP --parties $NB --logfile $TAG.log & echo $!" noback)

# Waiting for protocol to complete
echo "waiting for completion"
CHCK=$(./run_EC2.sh 1 "kill -0 $PID" 2>&1)

while [ "$CHCK" == "" ]; do
    CHCK=$(./run_EC2.sh 1 "kill -0 $PID" 2>&1)
    sleep 5
    done

./run_EC2.sh 1 "aws s3 cp /home/ubuntu/$TAG.log s3://ligero-simulations/" noback
