"""
Original Developer: David Roberts
Purpose of Module: solves time-dependent Redfield equation encoded
in the Linblads imported from linblad.csv
Last Modified: 6/5/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
from qutip import *
import numpy as np
import time

# Custom Modules:
import parameters
import helper
import ODE_integrator


globalstart = time.time()

# Load parameters
list_of_t = parameters.LIST_OF_TIMES
time_indices = range(len(list_of_t))
initial_density_matrix = helper.vectorize(parameters.INITIAL_DENSITY_MATRIX)

root_filepath = parameters.ROOT_FILEPATH
filename = parameters.FILENAME



print ("importing linbladians...")
start = time.time()
list_of_linblads_reals = np.loadtxt(root_filepath+filename+'reals_v2.csv', delimiter=",")
list_of_linblads_imags = np.loadtxt(root_filepath+filename+'imags_v2.csv', delimiter=",")
end = time.time()
print ("imported linbladians in {} seconds.".format(end-start))

# print(list_of_linblads_imags)

# print("list_of_linblads_imags[0] is {}".format(list_of_linblads_imags[0]))
# print("list_of_linblads_imags[1] is {}".format(list_of_linblads_imags[1]))


print ("combining real and imaginary parts...")
start = time.time()
list_of_linblads = np.array([linblad[0] + 1j * linblad[1] for linblad in 
									zip(list_of_linblads_reals, list_of_linblads_imags)])
end = time.time()
print ("combined real and imaginary parts in {} seconds.".format(end-start))

print ("reshaping linbladians...")
start = time.time()
dim = parameters.NUM_STATES_CUTOFF
linblad_operators = map(helper.get_compact_tensor_from_compact_compact_tensor,
										zip(list_of_linblads, [dim]*len(list_of_t)))
end = time.time()
print ("reshaped linbladians in {} seconds.".format(end-start))


print ("running time evolution...")
output_evolution = ODE_integrator.run_time_evolution(
	initial_density_matrix, linblad_operators, list_of_t, "TR")
# output_evolution = ODE_integrator.run_time_evolution(
# 	initial_density_matrix, linblad_operators, list_of_t, "Explicit Euler")

print ("completed time evolution in {} seconds. Exporting as csv...".format(end-start))

output_evolution_reals = np.array(map(np.real, output_evolution))
output_evolution_imags = np.array(map(np.imag, output_evolution))

np.savetxt(root_filepath + filename + "density_matrix_reals.csv", output_evolution_reals, delimiter=",")
np.savetxt(root_filepath + filename + "density_matrix_imags.csv", output_evolution_imags, delimiter=",")

np.savetxt(root_filepath + filename + "simulation_times.csv", list_of_t, delimiter=",")

globalend = time.time()
print ("td_linblad_solver.py complete. process took {} seconds.".format(globalend-globalstart))


