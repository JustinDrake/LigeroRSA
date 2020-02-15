#!/bin/bash

if [ $# -ne 1 ]; then
    echo "USAGE: propagate_keys <key_name>" 
    exit 1
fi

chmod 400 AWS/$1.pem
ssh-keygen -y -f AWS/$1.pem > AWS/$1.pub 

while IFS=, read -r region location_description instance_type nb_instances ami security_group
       do
 aws ec2 import-key-pair --key-name $1 --public-key-material file://$PWD/AWS/$1.pub --region $region
 done < AWS/map_deployment_rsa_ceremony
