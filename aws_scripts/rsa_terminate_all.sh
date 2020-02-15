#!/bin/bash

# To terminate all instances across all simulations, specify ALL, otherwise only instances for the current TAG will be terminated
if [ $# -gt 1 ]; then
    echo "USAGE: rsa_terminate_all.sh <optional:ALL>"
    exit 1
fi

# Apply filter if needed
if [ $# -lt 1 ]; then
    FILTER="Name=tag:role,Values=$TAG"
    echo "Terminating all instances with $FILTER"
else
    FILTER=""
    echo "Terminating all instances, no filter"
fi

# Proceed with terminating instances
while IFS=, read -r region location_description instance_type nb_instances ami security_group
    do 
        echo "Terminating Region $region"

        aws ec2 describe-instances --region $region --filters "Name=instance-state-name,Values=running" $FILTER | jq -r .Reservations[].Instances[].InstanceId | xargs -L 1 -I {} aws ec2 terminate-instances --region $region --instance-id {} &
    done < AWS/map_deployment_rsa_ceremony