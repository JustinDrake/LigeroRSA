#!/bin/bash

if [ $# -lt 1 ]; then

    NB=0
    while IFS=, read -r location location_description instance_type nb_instances ami security_group
        do    
                NB=$((NB+$nb_instances))
        done < AWS/map_deployment_rsa_ceremony
    echo "No number of parties specified, will use $NB parties based on configuration file"
 
else 
    NB=$1
fi

PUBLIC_IP=$(cat AWS/coordinator.ip)

instance_id=1
maxval=$((NB))

buffered_pids=()

while [ $instance_id -le $maxval ]; do
    LOCAL_IP=$(cat "AWS/Registered_IPs/$instance_id")
    # full protocol
    ./run_EC2.sh "$instance_id" "gunzip -f ./$PARTY_BINARY.gz ; ./$PARTY_BINARY --ip $PUBLIC_IP --localIPaddress $instance_id > log_stdout &" &

    buffered_pids+=($!)

    echo "$((instance_id-1)) participating."
    if (( $instance_id % $INDIVIDUAL_BATCH_SIZE == 0 ))
    then 
        # Wait for the deployment to be completed
        for buffered_pid in ${buffered_pids[*]}; do
            wait $buffered_pid
            done
        buffered_pids=()
    fi

    let "instance_id++"
    done

# Wait for the deployment to be completed
for buffered_pid in ${buffered_pids[*]}; do
    wait $buffered_pid
    done