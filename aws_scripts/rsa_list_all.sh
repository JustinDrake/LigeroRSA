#!/bin/bash

# Reading Command Line
POSITIONAL=()
ALERT=""
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -a|--alert)
    ALERT=true
    THRESHOLD="$2"
    shift                   # past argument
    shift                   # past value
    ;;
    *)                      # unknown option
    POSITIONAL+=("$1")      # save it in an array for later
    shift                   # past argument
esac
done
set -- "${POSITIONAL[@]}"   # restore positional parameters
    
if [ $# -gt 0 ]; then
    echo "USAGE: rsa_list_all.sh"
    echo
    echo "-a, --alert <threshold>: sends an email alert when the number of instances is above the threshold"
    exit 1
fi

# Initialize parameters
if [ -z "$TAG" ]; then
    FILTER=""
    echo "No tag filter."
else
    FILTER="Name=tag:role,Values=$TAG"
    echo "Filtering instances with tag $FILTER"
fi

echo $FILTER

# Preparing Report
echo "AWS Usage Report By Region:" > AWS/usage.csv
echo "location, EC2 region code, # instances" >> AWS/usage.csv

# Checking EC2 Regions
echo
echo "Querying EC2 Regions"
echo "--------------------"

while IFS=, read -r region location_description instance_type nb_instances ami security_group
do                
    # Interprocess communication pipe
    PIPE_RESULTS="/tmp/aws_usage_list_$$_$region"
    mkfifo $PIPE_RESULTS

    (
        NB=$(aws ec2 describe-instances --filters "Name=instance-state-name,Values=running" $FILTER --region $region --query 'Reservations[*].Instances[*].[PublicDnsName]' --output text |  wc -l) 
        echo "$location_description,$region,$NB" > $PIPE_RESULTS
    ) &
done < AWS/map_deployment_rsa_ceremony

IDX=0
NB=0

while IFS=, read -r region location_description instance_type nb_instances ami security_group
do                
    PIPE_RESULTS="/tmp/aws_usage_list_$$_$region"
    LINE=$(< $PIPE_RESULTS)
    echo $LINE >> AWS/usage.csv
    NB_ITER=$(echo $LINE | tr "," "\n" | tail -n1)
    NB=$(( $NB_ITER + $NB ))
done < AWS/map_deployment_rsa_ceremony

# Finalizing Report
echo "Aggregate,aggregate,$NB" >> AWS/usage.csv

# Display the Report
cat AWS/usage.csv | column -t -s,

# Send alert if required
if [ "$ALERT" = true ] && [ $NB -ge $THRESHOLD ]; then
    echo "alert! sending email"
    ./ligero_mail_alert.sh
fi
