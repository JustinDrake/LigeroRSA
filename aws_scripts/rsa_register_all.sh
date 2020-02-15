#!/bin/bash

if [ $# -ne 0 ]; then
    echo "USAGE: rsa_register_all.sh"
    exit 1
fi

# Clean up existing records
rm -rf AWS/Registered_IPs/*
rm -rf AWS/Registered_Private_IPs/*
rm -rf AWS/Registered_Instances/*

# Register all
while IFS=, read -r location_code location_description instance_type nb_instances ami security_group
do        
    echo "Registering $location_code"
    ./rsa_register_region.sh $location_code
done < AWS/map_deployment_rsa_ceremony
echo "AWS EC2 Provisioning Completed."