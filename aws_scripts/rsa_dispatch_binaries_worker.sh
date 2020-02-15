#!/bin/bash

# Establishing the number of parties
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

# Clean up
rm -rf logs/*

# Record of all detached processes
process_pids=()
buffered_pids=()

PUBLIC_IP=$(cat AWS/Registered_IPs/1)
PRIVATE_IP=$(cat AWS/Registered_Private_IPs/1)

instance_id=1
maxval=$(($NB))

# Compress the Party Binary
gzip -fk "$BUILD_DIR/$PARTY_BINARY"

while [ $instance_id -le $maxval ]; do

    # Dispatching the parties
    ./put_EC2.sh "$instance_id" "$BUILD_DIR/$PARTY_BINARY.gz" / 2> logs/"$instance_id"_coordinator_error_log & 
    process_pids+=($!)  
    buffered_pids+=($!)  

    ./put_EC2.sh "$instance_id" AWS/logging / 2> logs/"$instance_id"_logging_error_log & 
    process_pids+=($!) 
    buffered_pids+=($!)  

    echo "Dispatching party $((instance_id-1))" 
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

# Waiting for full deployment to complete
echo "Waiting for deployment completion."

# Wait for the deployment to be completed
for process_pid in ${process_pids[*]}; do
    wait $process_pid
    done

# Clean up 
rm -f "$BUILD_DIR/$PARTY_BINARY.gz"

echo "Deployment complete."