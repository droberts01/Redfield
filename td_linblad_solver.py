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


# Custom Modules:
import parameters
import helper

# Load parameters
list_of_t = parameters.LIST_OF_TIMES
time_indices = range(len(list_of_t))
initial_density_matrix = helper.vectorize(parameters.INITIAL_DENSITY_MATRIX)


list_of_linblads_reals = np.loadtxt('data/linblad_real_v2.csv', delimiter=",")
list_of_linblads_imags = np.loadtxt('data/linblad_imag_v2.csv', delimiter=",")


list_of_linblads = np.array([linblad[0] + 1j * linblad[1] for linblad in 
									zip(list_of_linblads_reals, list_of_linblads_imags)])


dim = parameters.NUM_STATES_CUTOFF
linblad_operators = map(helper.get_compact_tensor_from_compact_compact_tensor,
										zip(list_of_linblads, [dim]*len(list_of_t)))


output_evolution = [initial_density_matrix]*len(time_indices)

for time_index in time_indices:
	if time_index == 0:
		pass
	else:
		time_step = list_of_t[time_index] - list_of_t[time_index - 1]
		if time_index < 150:
			print("time_index is {}".format(time_index))
			print("output_evolution[time_index-1] is {}".format(output_evolution[time_index-1]))
		incremental_evolution = time_step * np.matmul(linblad_operators[time_index - 1], 
												output_evolution[time_index - 1])
		output_evolution[time_index] = output_evolution[time_index - 1] + incremental_evolution 
			


output_evolution_reals = np.array(map(np.real, output_evolution))
output_evolution_imags = np.array(map(np.imag, output_evolution))

np.savetxt("/Users/Droberts/Documents/LANLA/Redfield/data/density_matrix_reals.csv", output_evolution_reals, delimiter=",")
np.savetxt("/Users/Droberts/Documents/LANLA/Redfield/data/density_matrix_imags.csv", output_evolution_imags, delimiter=",")

np.savetxt("/Users/Droberts/Documents/LANLA/Redfield/data/simulation_times.csv", list_of_t, delimiter=",")




