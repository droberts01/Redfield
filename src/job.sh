#!/bin/bash
# Shell script for automating job schedules


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
for N in 9
do
	for J in .3
	do 
		python generate.py 5E-6 .2 $J 1 $N 4 0.001 10 10000 1 "Home"


# solve.py -tQA, -I, -J, -K, -N, -Nc, 
#   			-step, -window_size, -num_samples, -decoherence
# Solve the Bloch-Redfield master equation (via 2nd order Implicit Runge-Kutta):
		python solve.py 5E-6 .2 $J 1 $N 4 0.001 10 10000 1 "Home"
	done
done

# for D in 0
# do

# done


ENDTIME=$(date +%s)




echo "Stopping:"
date

echo "Redfield simulation ran in $(($ENDTIME - $STARTTIME)) seconds."