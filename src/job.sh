#!/bin/bash
# My first script

echo "Checking:"
pwd
hostname
date
env | sort > ENVIRONMENT

echo "Starting Redfield simulation of the D Wave 2X at LANL:" 

# module load anaconda

STARTTIME=$(date +%s)
# generate.py -tQA, -I, -J, -K, -N, -Nc, 
#   			-step, -window_size, -num_samples, -decoherence
python generate.py 5E-6 .2 .23 1 5 4 10 0.001 20 10000 0

# for D in 0
# do
# python linblad.py 5 .2 .23 4 $D
# done

ENDTIME=$(date +%s)


echo "Stopping:"
date

echo "Redfield simulation ran in $(($ENDTIME - $STARTTIME)) seconds."