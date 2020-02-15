#!/bin/bash

if [ $# -ne 1 ]; then
    echo "USAGE: ec2_rsa_ceremony_deploy.sh <d(deploy)/p(provision)>"
    exit 1
fi

if [ -f "AWS/$EC2_KEY_ID.pem" ]; then

    if [ "$1" == "d" ]; then
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

    else 
        # Populating new records
        ./rsa_register_all.sh

        PUBLIC_IP=$(cat AWS/Registered_IPs/1)

        # Dispatching the coordinator 
        echo "Dispatching coordinator"
        ./put_EC2.sh 1 "$BUILD_DIR/$COORDINATOR_BINARY" /
        ./put_EC2.sh 1 "AWS/logging" /
        ./put_EC2.sh 1 "AWS/config" /

        # echo "Provisioning maps"
        # Instructions for each worker
        echo 2
        ./put_EC2.sh 2 AWS/map_deployment_rsa_ceremony.2 Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/AWS/map_deployment_rsa_ceremony
        echo 3
        ./put_EC2.sh 3 AWS/map_deployment_rsa_ceremony.3 Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/AWS/map_deployment_rsa_ceremony
        echo 4
        ./put_EC2.sh 4 AWS/map_deployment_rsa_ceremony.4 Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/AWS/map_deployment_rsa_ceremony
        echo 5
        ./put_EC2.sh 5 AWS/map_deployment_rsa_ceremony.5 Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/AWS/map_deployment_rsa_ceremony

        # Send over coordinator's public IP
        echo "Passing on coordinator ip"
        echo $PUBLIC_IP > "coordinator.ip"
        echo 2
        ./put_EC2.sh 2 coordinator.ip Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/AWS/
        echo 3
        ./put_EC2.sh 3 coordinator.ip Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/AWS/
        echo 4
        ./put_EC2.sh 4 coordinator.ip Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/AWS/
        echo 5
        ./put_EC2.sh 5 coordinator.ip Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/AWS/

        # Copy binaries
        echo "Provisioning binaries"
        echo 2
        ./put_EC2.sh 2 ../build/src/party_full_protocol Documents/Ligero/code/current/yuriy/RSACeremony/build/src/
        ./put_EC2.sh 2 ../build/src/coordinator_full_protocol Documents/Ligero/code/current/yuriy/RSACeremony/build/src/
        echo 3
        ./put_EC2.sh 3 ../build/src/party_full_protocol Documents/Ligero/code/current/yuriy/RSACeremony/build/src/
        ./put_EC2.sh 3 ../build/src/coordinator_full_protocol Documents/Ligero/code/current/yuriy/RSACeremony/build/src/
        echo 4
        ./put_EC2.sh 4 ../build/src/party_full_protocol Documents/Ligero/code/current/yuriy/RSACeremony/build/src/
        ./put_EC2.sh 4 ../build/src/coordinator_full_protocol Documents/Ligero/code/current/yuriy/RSACeremony/build/src/
        echo 5
        ./put_EC2.sh 5 ../build/src/party_full_protocol Documents/Ligero/code/current/yuriy/RSACeremony/build/src/
        ./put_EC2.sh 5 ../build/src/coordinator_full_protocol Documents/Ligero/code/current/yuriy/RSACeremony/build/src/

        # Send workers scripts
        echo "Provisioning infrastructure scripts"
        echo 2
        ./put_EC2.sh 2 rsa_ceremony_deploy_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        ./put_EC2.sh 2 rsa_ceremony_start_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        ./put_EC2.sh 2 rsa_dispatch_binaries_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        echo 3
        ./put_EC2.sh 3 rsa_ceremony_deploy_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        ./put_EC2.sh 3 rsa_ceremony_start_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        ./put_EC2.sh 3 rsa_dispatch_binaries_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        echo 4
        ./put_EC2.sh 4 rsa_ceremony_deploy_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        ./put_EC2.sh 4 rsa_ceremony_start_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        ./put_EC2.sh 4 rsa_dispatch_binaries_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        echo 5
        ./put_EC2.sh 5 rsa_ceremony_deploy_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        ./put_EC2.sh 5 rsa_ceremony_start_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
        ./put_EC2.sh 5 rsa_dispatch_binaries_worker.sh Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts/
    fi

else 
    echo "Error: AWS private key missing, aborting now."
fi