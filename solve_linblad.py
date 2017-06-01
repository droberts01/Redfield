"""
Original Developer: David Roberts
Purpose of Module: outputs solution to a time-dependent Linblad equation.
Last Modified: 5/30/17
Last Modified By: David Roberts
Last Modification Purpose: fixed function naming
"""


# Standard Modules:
from qutip import *
import csv
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as linalg
import pandas as pd


# Custom Modules:
import parameters
import redfield
import primitive
import linblad_solver


# Load global variables
TIME_STEP = parameters.TIME_STEP
NUM_STATES_CUTOFF = parameters.NUM_STATES_CUTOFF
ANNEALING_PARAMETER = parameters.ANNEALING_PARAMETER
INITIAL_DENSITY_MATRIX = parameters.INITIAL_DENSITY_MATRIX

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def tensor_apply(linbladian, density_matrix, s):

	def output_components_before_contraction(i,j,k,l):
		return linbladian[s][i][j][k][l]*density_matrix[i][j]

	def output_components(k,l):
		output = 0
		for i in range(NUM_STATES_CUTOFF):
			for j in range(NUM_STATES_CUTOFF):
				output = output + output_components_before_contraction(i,j,k,l)
		return output

	return [[output_components(k,l) for l in range(NUM_STATES_CUTOFF)] for k in range(NUM_STATES_CUTOFF)]


def compute_linblad_evolution(initial_density_matrix, linbladian):
	density_matrix_linblad_evolution = (([0]*NUM_STATES_CUTOFF)*NUM_STATES_CUTOFF)*len(ANNEALING_PARAMETER)
	transient_density_matrix = initial_density_matrix
	density_matrix_linblad_evolution[0] = initial_density_matrix
	for s in ANNEALING_PARAMETER:
		infinitessimal_change = tensor_apply(linbladian, transient_density_matrix, s)
		transient_density_matrix = transient_density_matrix + TIME_STEP*infinitessimal_change
		density_matrix_linblad_evolution[s] = transient_density_matrix
	return density_matrix_linblad_evolution



# Compute linbladian
linblad_tensor = linbladian
time_dependent_density_matrix = compute_linblad_evolution(INITIAL_DENSITY_MATRIX)
time_dependent_density_matrix_array_version = [sum([time_dependent_density_matrix[s][i] for i in range(NUM_STATES_CUTOFF)]) for s in ANNEALING_PARAMETER]

dict_list = [{'s'+str(s) : time_dependent_density_matrix_array_version[s]} for s in range(len(ANNEALING_PARAMETER]))]
column_list = ['s'+str(s) for s in range(len(ANNEALING_PARAMETER))]

final_dict = {}
for dictionary in dict_list:
	final_dict = merge_two_dicts(final_dict, dictionary)


df = pd.DataFrame(final_dict, columns = column_list)



