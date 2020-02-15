#!/bin/bash

if [ $# -ne 0 ]; then
    echo "USAGE: ec2_rsa_ceremony_deploy.sh"
    exit 1
fi

if [ -f "AWS/$EC2_KEY_ID.pem" ]; then

    NBTOTALINSTANCES=0
    while IFS=, read -r location_code location_description instance_type nb_instances ami security_group
        do
            echo "Spawning $nb_instances instance(s) of AMI $ami in $location_description (AWS $location_code), security group = $security_group."
            NBTOTALINSTANCES=$((NBTOTALINSTANCES + $nb_instances))
            aws ec2 run-instances --no-dry-run --region=$location_code --image-id "$ami" --count $nb_instances --instance-type "$instance_type" --key-name "$EC2_KEY_ID" --tag-specifications "ResourceType=instance,Tags=[{Key=role,Value=$TAG}]" --security-group-ids "$security_group" > /dev/null &
        done < AWS/map_deployment_rsa_ceremony

    # Wait until all instances are ready
    ACTUAL=0
    while [ $ACTUAL -lt $(( NBTOTALINSTANCES )) ]; do
        NB=0
        while IFS=, read -r location location_description instance_type nb_instances ami security_group
            do
            NB1=$(aws ec2 describe-instances --region $location --filters "Name=instance-state-name,Values=running" "Name=tag:role,Values=$TAG" --query 'Reservations[*].Instances[*].[PublicDnsName]' --output text | wc -l)
            NB=$((NB+NB1))

            done < AWS/map_deployment_rsa_ceremony
            ACTUAL=$NB
        
        if [ $ACTUAL -lt $(( NBTOTALINSTANCES )) ]; then 
            echo "Only $ACTUAL instances out of $NBTOTALINSTANCES ready, still sleeping."
            sleep $EC2_STATUS_CHECK_FREQUENCY
        else
            echo "Deployment complete."
        fi
    done

    # Populating new records
    ./rsa_register_all.sh

    # Dispatch binaries
    ./rsa_dispatch_binaries_worker.sh $((NBTOTALINSTANCES))

    echo "Configuration Ready for RSA Ceremony."

else 
    echo "Error: AWS private key missing, aborting now."
fi