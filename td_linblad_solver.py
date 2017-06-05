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
import numpy
import pandas as pd


# Custom Modules:
import parameters
import helper

# Load parameters
list_of_t = parameters.LIST_OF_TIMES
time_indices = range(len(list_of_t))
initial_density_matrix = parameters.INITIAL_DENSITY_MATRIX


list_of_linblads_csv = pd.read_csv('linblad.csv', sep=',',header=None)
list_of_linblads = map(helper.init_compact_compact_tensor_from_array, list_of_linblads_csv)
list_of_linblad_operators = map(helper.get_compact_tensor_from_compact_compact_tensor, list_of_linblads)
linblad_operators = [linblad.array for linblad in list_of_linblad_operators]


output_evolution = initial_density_matrix*len(time_indices)
for time_index in time_indices:
	if time_index == 0:
		pass
	else:
		time_step = list_of_t[time_index] - list_of_t[time_index - 1]
		incremental_evolution = time_step * np.matmul(linblad_operators[time_index - 1], 
												output_evolution[time_index - 1])
		output_evolution[time_index] = output_evolution[time_index - 1] + incremental_evolution 

output_evolution = np.array(output_evolution)
np.savetxt("redfield_simulation.csv", output_evolution, delimiter=",")
np.savetxt("redfield_simulation_times.csv", list_of_t, delimiter=",")




