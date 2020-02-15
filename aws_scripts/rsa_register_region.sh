#!/bin/bash

if [ $# -ne 1 ]; then
    echo "USAGE: rsa_register_region.sh <region>"
    exit 1
fi

# Gather all current EC2 instances
FILTER="Name=tag:role,Values=$TAG"
PUBLIC_IP=$(aws ec2 describe-instances --region $1 --filters "Name=instance-state-name,Values=running" $FILTER --query 'Reservations[*].Instances[*].[PublicIpAddress]' --output text)
PRIVATE_IP=$(aws ec2 describe-instances --region $1 --filters "Name=instance-state-name,Values=running" $FILTER --query 'Reservations[*].Instances[*].[PrivateIpAddress]' --output text)
DNS=$(aws ec2 describe-instances --region $1 --filters "Name=instance-state-name,Values=running" $FILTER --query 'Reservations[*].Instances[*].[PublicDnsName]' --output text)

# Populate the IPs
REF=$(ls AWS/Registered_IPs/ -1 | wc -l)

for i in $(echo $PUBLIC_IP)
do
REF=$((REF+1))
echo $i > "AWS/Registered_IPs/$REF"
done

REF=$(ls AWS/Registered_Private_IPs/ -1 | wc -l)

for i in $(echo $PRIVATE_IP)
do
REF=$((REF+1))
echo $i > "AWS/Registered_Private_IPs/$REF"
done

# Populate the DNSs
REF=$(ls AWS/Registered_Instances/ -1 | wc -l)

for i in $(echo $DNS)
do
REF=$((REF+1))
echo "ubuntu@$i" > "AWS/Registered_Instances/$REF"
done
