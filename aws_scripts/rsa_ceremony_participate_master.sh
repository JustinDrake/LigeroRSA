#!/bin/bash

if [ $# -ne 0 ]; then
    echo "USAGE: ec2_rsa_ceremony_deploy.sh"
    exit 1
fi

if [ -f "AWS/$EC2_KEY_ID.pem" ]; then

    echo "Executing"
    ./run_EC2.sh 2 "cd ~/Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts;source AWS/config Tarik w_2;./rsa_ceremony_deploy_worker.sh;./rsa_ceremony_start_worker.sh" > log_worker_2 &
    ./run_EC2.sh 3 "cd ~/Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts;source AWS/config Tarik w_3;./rsa_ceremony_deploy_worker.sh;./rsa_ceremony_start_worker.sh" > log_worker_3 &
    ./run_EC2.sh 4 "cd ~/Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts;source AWS/config Tarik w_4;./rsa_ceremony_deploy_worker.sh;./rsa_ceremony_start_worker.sh" > log_worker_4 &
    ./run_EC2.sh 5 "cd ~/Documents/Ligero/code/current/yuriy/RSACeremony/aws_scripts;source AWS/config Tarik w_5;./rsa_ceremony_deploy_worker.sh;./rsa_ceremony_start_worker.sh" > log_worker_5 &

else 
    echo "Error: AWS private key missing, aborting now."
fi