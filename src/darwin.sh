#!/bin/bash
# Shell script for automating job schedules
#SBATCH --job-name=comp1
#SBATCH --mail-user=youremail@gmail.com
#SBATCH --mail-type=END
#SBATCH --nodes=1 
#SBATCH --exclusive 
#SBATCH --time=1-00:00:00 
#SBATCH --no-requeue 
#SBATCH --output=out 
#SBATCH --error=err 
#SBATCH --mem 20000


set -x                          # Output commands 
set -e                          # Abort on errors 



echo "Checking:"
pwd
hostname
date
env | sort > ENVIRONMENT

echo "Starting Redfield simulation of the D Wave 2X:" 


# Needed on certain supercomputing clusers:
# module load anaconda



STARTTIME=$(date +%s)


# generate.py -tQA, -I, -J, -K, -N, -Nc, 
#   			-step, -window_size, -num_samples, -decoherence
# Generate a Bloch-Redfield master equation:

python generate.py 5E-6 .2 .24 1 7 4 0.001 10 1200 0 "Darwin"


# solve.py -tQA, -I, -J, -K, -N, -Nc, 
#   			-step, -window_size, -num_samples, -decoherence
# Solve the Bloch-Redfield master equation (via 2nd order Implicit Runge-Kutta):
python solve.py 5E-6 .2 .24 1 7 4 0.001 10 1200 0 "Darwin"



# for D in 0
# do

# done


ENDTIME=$(date +%s)




echo "Stopping:"
date

echo "Redfield simulation ran in $(($ENDTIME - $STARTTIME)) seconds."