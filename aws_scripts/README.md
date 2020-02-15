# AWS Scripts for RSA Ceremony 

## Run a full simulation (with automated provisioning, executing and results exported to S3)
Commands are accessible from the aws_script folder:
- native_build.sh: Build the binary for the first time (including fetches dependencies)
- rebuild_and_run_simulation.sh: Rebuild the binaries & run simulation
- run_simulation.sh: Just run the simulation

## Execute a dry-run
- Pre-requisite: source AWS/config any_tag
- This can be helpful when you want to check that you have the proper credentials, enough instances available on EC2 etc. to run the simulation, without actually consuming any resources. It will attempt to provision on EC2 based on current configuration.
- If successful, the following error will be produced: An error occurred (DryRunOperation) when calling the RunInstances operation: Request would have succeeded, but DryRun flag is set.

## Run a simulation (manually, step by step)
Source AWS/config (source AWS/config <tag>) to set up environment variables beforehand, then run:
   ./rsa_ceremony_deploy.sh
   ./rsa_ceremony_start.sh
   ./rsa_terminate_all.sh 

## View, terminate running instances
Source AWS/config (source AWS/config <tag>) to set up environment variables beforehand, then run:
   ./rsa_list_all.sh: view all running instances with tag <tag>
   ./rsa_terminate_all.sh: terminate all instances with tag <tag>

   - Terminating a single instance can be done via the AWS EC2 web portal.
   - Set TAG to any other value to switch to another live simulation.
   - Set TAG="" to remove the filter by simulation tag, and view/terminate ALL EC2 instances (use with caution).

## Configuration files

### AWS/config:
    This contains all parameters for the simulation as well as configuration options such as the EC2 key ID being used. 
    
### AWS/map_deployment_rsa_ceremony:
    The configuration parameters for the simulation are stored in AWS/map_deployment_rsa_ceremony. The format of the file is:
    <EC2 location code> <region full description> <instance type> <nb of instances> <Amazon AMI deployed> <Security group>

### AWS/mail_alerts:
    List of email addresses for alerts. 

### AWS/logging:
    All parameters for logging during the simulation

## Scripts for RSA Manual Simulation

### rsa_ceremony_start.sh:
    Kicks off the ceremony. Parties start up first (for many parties, they will start in batches - the size of the individual batches can be configured in this script in the parameter INDIVIDUAL_BATCH_SIZE)

### rsa_ceremony_deploy.sh:
    This script spawns instances as specified in the configuration file, and once they are active registers them for future use (notes down IP and instance name)

### rsa_ceremony_deploy_dry_run.sh:
    Use this to test AWS EC2 CLI authorization for a configuration file. No instance will be spawned and no infrastructure deployed. If successful, the test will return error messages indicating that the functions could not be run because we are in dry_run mode.

### rsa_register_all.sh:
    Registers IPs and instance names for all active instances across all regions described in the configuration file.

### rsa_register_region.sh:
    Registers IPs and instance names for all instances within a specific region.

### rsa_list_all.sh:
    Lists all active instances across all regions included in the configuration file, sends an email alert if the instance count is above the threshold.

### rsa_dispatch_binaries.sh:
    Dispatch the actual binaries for the parties and coordinator to the appropriate EC2 instances. EC2 instance 1 hosts the coordinator, EC2 instances 2 to N host the parties

### rsa_terminate_all.sh:
    Identifies all active instances across all the regions included in the configuration file, and terminates them.

## Administative scripts

### ligero_mail_alert.sh:
    Sends an e-mail alert.

### go_EC2.sh:
    Establishes an SSH connection to an EC2 instance using its internal reference number.

### run_EC2.sh:
    Runs a command on an EC2 instance based on its internal reference number.

### put_EC2.sh:
    Pushes out a file to an EC2 instance using its internal reference number.

### propagate_keys.sh:
    Generates a public key based on the private key file specified and exports it to all AWS regions.
