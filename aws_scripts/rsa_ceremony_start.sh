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

PUBLIC_IP=$(cat AWS/Registered_IPs/1)
PRIVATE_IP=$(cat AWS/Registered_Private_IPs/1)

instance_id=2
maxval=$((NB+1))

buffered_pids=()

while [ $instance_id -le $maxval ]; do
    LOCAL_IP=$(cat "AWS/Registered_IPs/$instance_id")
    # full protocol
    ./run_EC2.sh "$instance_id" "gunzip -f ./$PARTY_BINARY.gz ; ./$PARTY_BINARY --ip $PUBLIC_IP --localIPaddress $LOCAL_IP > log_stdout &" > /dev/null &
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

# Kick off protocol
./put_EC2.sh 1 "AWS/config" /
# full protocol
PID=$(./run_EC2.sh 1 "(source config $TAG) && (ulimit -n $COORDINATOR_NB_CONCURRENT_CONNECTIONS) && (./$COORDINATOR_BINARY --logfile $TAG.log --parties $NB --ip $PRIVATE_IP > log_coordinator) & echo $!" noback 2> /dev/null)

# Waiting for protocol to complete
echo "waiting for completion"
CHCK=$(./run_EC2.sh 1 "kill -0 $PID" 2>&1)

while [ "$CHCK" == "" ]; do
    CHCK=$(./run_EC2.sh 1 "kill -0 $PID" 2>&1)
    sleep 5
    done

./run_EC2.sh 1 "aws s3 cp /home/ubuntu/$TAG.log s3://ligero-simulations/" noback
