#!/bin/bash

if [ $# -ne 0 ]; then
    echo "USAGE: ec2_rsa_ceremony_deploy_dry_run.sh"
    exit 1
fi

NBTOTALINSTANCES=0
while IFS=, read -r location_code location_description instance_type nb_instances ami security_group
    do    
        echo "Spawning $nb_instances instance(s) of AMI $ami in $location_description (AWS $location_code), security group = $security_group."
        NBTOTALINSTANCES=$((NBTOTALINSTANCES + $nb_instances))
    aws ec2 run-instances --dry-run --tag-specifications 'ResourceType=instance,Tags=[{Key=role,Value=$TAG}]' --region=$location_code --image-id "$ami" --count $nb_instances --instance-type "$instance_type" --key-name "$EC2_KEY_ID" --security-group-ids "$security_group" > /dev/null &
    done < AWS/map_deployment_rsa_ceremony

echo "Total # of instances: $NBTOTALINSTANCES"
