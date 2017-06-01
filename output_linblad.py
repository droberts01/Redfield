"""
Original Developer: David Roberts
Purpose of Module: outputs Redfield Linblad.
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


# Load global variables
MODE = parameters.MODE
NUM_STATES = parameters.NUM_STATES
NUM_STATES_CUTOFF = parameters.NUM_STATES_CUTOFF
TIME_STEP = parameters.TIME_STEP
ANNEALING_PARAMETER = parameters.ANNEALING_PARAMETER


def diabatic_tensor(eigenvectors):
    states = range(NUM_STATES_CUTOFF)
    def basis_dragging_matrix_components(s,i,j):
        eigenbra_derivative = eigenvectors[s + 1][i] + eigenvectors[s - 1][i]
        eigenket_derivative = eigenvectors[s + 1][i] - eigenvectors[s - 1][i]
        identity_matrix = identity(NUM_STATES)
        braket_derivative = identity_matrix.matrix_element(eigenbra_derivative, eigenket_derivative)/(4*TIME_STEP)
        return braket_derivative
    
    def diabatic_tensor_components(s,i,j,k,l):
        if s == 0 or s == len(ANNEALING_PARAMETER)-1:
            return 0
        else:
            return helper.delta(i,k)*basis_dragging_matrix_components(s,l,j) + helper.delta(j,l)*basis_dragging_matrix_components(s,i,k)
    
    return [[[[[diabatic_tensor_components(s,i,j,k,l) for l in range(NUM_STATES_CUTOFF)] for k in range(NUM_STATES_CUTOFF)] for j in range(NUM_STATES_CUTOFF)] for i in range(NUM_STATES_CUTOFF)] for s in ANNEALING_PARAMETER]


def frequency_tensor(frequency_matrix):
    def frequency_tensor_components(s,i,j,k,l):
        return helper.delta(j,l)*helper.delta(i,k)*frequency_matrix[s][i][j]
    return [[[[[frequency_tensor_components(s,i,j,k,l) for l in range(NUM_STATES_CUTOFF)] for k in range(NUM_STATES_CUTOFF)] for j in range(NUM_STATES_CUTOFF)] for i in range(NUM_STATES_CUTOFF)] for s in ANNEALING_PARAMETER]



def linbladian:

	# Initialize and diagonalize D Wave Hamiltonians
    if MODE == "GOOGLE":
        hamiltonian = [google_quail_primitive.d_wave_hamiltonian(s) for s in ANNEALING_PARAMETER]
    elif MODE == "LANL":
        hamiltonian = [primitive.d_wave_hamiltonian(s) for s in ANNEALING_PARAMETER]
    else:
        print "ERR"

    eigenvalues = [hamiltonian[s].eigenenergies() for s in ANNEALING_PARAMETER]
    eigenvectors = [hamiltonian[s].eigenstates() for s in ANNEALING_PARAMETER]
    idx = [eigenvalues[s].argsort()[::-1] for s in ANNEALING_PARAMETER]
    sorted_eigenvalues = [eigenvalues[s][idx[s]] for s in ANNEALING_PARAMETER]
    sorted_eigenvectors = [eigenvectors[s][:,idx[s]] for s in ANNEALING_PARAMETER]
    truncated_eigenvalues = eigenvalues_sorted[::NUM_STATES_CUTOFF]
    
    # Raw spectral data for computation of Redfield and diabatic tensors
    truncated_eigenvectors = sorted_eigenvectors[::STATE_CUT_OFF]
    frequency_matrix = [[sorted_eigenvalues[i]-sorted_eigenvalues[j] for i in range(NUM_STATES)] for j in range(NUM_STATES)]
    truncated_frequency_matrix = frequency_matrix[::NUM_STATES_CUTOFF][::NUM_STATES_CUTOFF]
    
    # non-interacting term
    noninteracting_term = -1j*frequency_tensor(truncated_frequency_matrix)

    # Redfield term
    redfield_term = redfield.redfield_tensor(sorted_eigenvectors, frequency_matrix)
 
    # diabatic term
    diabatic_term = diabatic_tensor(truncated_eigenvectors)
    
    return noninteracting_term - (redfield_term - diabatic_term)
    

# 
