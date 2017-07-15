#!/bin/bash
# Shell script for automating job schedules
#----------------

echo "Checking:"
pwd
hostname
date
env | sort > ENVIRONMENT

echo "Starting Redfield simulation of the D Wave 2X:" 


# Needed on certain supercomputing clusers:
# module load anaconda


#----------------------------------START-----------------------~~~~~
STARTTIME=$(date +%s)

#-----------------------------------------------------------------
# generate.py -tQA, -I, -J, -K, -N, -Nc, 
#   			-step, -window_size, -num_samples, -decoherence
# Generate a Bloch-Redfield master equation:

python generate.py 5E-6 .2 .23 1 5 4 10 0.001 20 15000 0

#-----------------------------------------------------------------
# solve.py -tQA, -I, -J, -K, -N, -Nc, 
#   			-step, -window_size, -num_samples, -decoherence
# Solve the Bloch-Redfield master equation (via 2nd order Implicit Runge-Kutta):
python solve.py 5E-6 .2 .23 1 5 4 10 0.001 20 15000 0


#----------------------------------FOR LOOP-----------------------
# for D in 0
# do

# done


ENDTIME=$(date +%s)
#----------------------------------END-----------------------~~~~~



echo "Stopping:"
date

echo "Redfield simulation ran in $(($ENDTIME - $STARTTIME)) seconds."