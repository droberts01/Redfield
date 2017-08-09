#!/bin/bash

# module load anaconda

python generate.py 5E-6 .2 0.43 1 6 4 0.001 10 10000 1 1 0 'Home'
python solve.py 5E-6 .2 0.43 1 6 4 0.001 10 10000 1 1 0 'Home'
