#---------------
# solve.py 
# Solves the Redfield quantum master equation
from sys import argv
import json
import meta_functions
import integrate

# Determine which JSON file to import
terminal_input = argv[1:]
args = [float(num) for num in terminal_input]



# solve Redfield master equation
from time import time

start = time()
integrate.solve_json(args)
end = time()

print("Overall simulation runtime was {} seconds.".format(end-start))


# initial_condition = np.zeros((Nc, Nc))
# initial_condition[0,0] = 1


