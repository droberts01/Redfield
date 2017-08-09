# -*- coding: utf-8 -*-
import re
import os
import subprocess
from subprocess import PIPE
import time
import numpy as np

COUNTER = 0
Nvals = [6, 8, 10]
Jvals = np.arange(.24, .5, .04)
# Jvals = [.39, .41, .43]
PAUSE_TIME = 2


for N in Nvals:
	for J in Jvals:
		sh_filename = 'darwin'+str(COUNTER)+'.sbatch'
		with open(sh_filename, "w") as file:

			file.write("#!/bin/bash\n")
			file.write("# Shell script for automating job schedules\n")
			file.write("#SBATCH --job-name=comp1\n")
			file.write("#SBATCH --mail-user=david.b.roberts@outlook.com\n")
			file.write("#SBATCH --mail-type=END\n")
			file.write("#SBATCH --nodes=1\n")
			file.write("#SBATCH --exclusive\n")
			file.write("#SBATCH --time=1-00:00:00\n") 
			file.write("#SBATCH --no-requeue \n")
			file.write("#SBATCH --output=out \n")
			file.write("#SBATCH --error=err \n")
			file.write("#SBATCH --mem 50000\n")
			file.write("set -x                          # Output commands\n")
			file.write("set -e                          # Abort on errors\n")



			file.write("echo 'Checking:'\n")
			file.write("pwd\n")
			file.write("hostname\n")
			file.write("date\n")
			file.write("env | sort > ENVIRONMENT\n")

			file.write("echo 'Starting Redfield simulation of the D Wave 2X:'\n")


			file.write("module load anaconda\n")
			file.write('\n')
			file.write("python generate.py 5E-6 .2 "+str(J)+" 1 "+str(N)+" 12 0.0002 10 10000 1 1 0 'Darwin'\n")
			file.write("python solve.py 5E-6 .2 "+str(J)+" 1 "+str(N)+" 12 0.0002 10 10000 1 1 0 'Darwin'\n")
		
			file.write("ENDTIME=$(date +%s)")




			file.write("echo 'Stopping:'\n")
			file.write("date\n")

			file.write("echo 'Redfield simulation ran in $(($ENDTIME - $STARTTIME)) seconds.'\n")
		
		print(" writing shell script with settings (5E-6 .2 "+str(J)+" 1 "+str(N)+" 12 0.0002 10 10000 1 1 0 'Darwin')")
		os.system('chmod -R 777 '+ sh_filename)


		command = 'sbatch ' + sh_filename
		print("executing {}".format(command))
		os.system(command)


		time.sleep(PAUSE_TIME)
		COUNTER += 1
		print('Finished submitting job. Starting a new job...')


print("done.")
