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

# Custom Modules:
import parameters
import redfield
import primitive
import linblad_solver


# Load global variables
TIME_STEP = parameters.TIME_STEP
NUM_STATES_CUTOFF = parameters.NUM_STATES_CUTOFF
ANNEALING_PARAMETER = parameters.ANNEALING_PARAMETER




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








