"""
Original Developer: David Roberts
Purpose of Module: outputs redfield Linblad.
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


def diabatic_tensor(eigenvectors):
    states = range(len(eigenvectors[0]))
    def basis_dragging_matrix_components(s,i,j):
        eigenbra_derivative = eigenvectors[s + 1][i] + eigenvectors[s - 1][i]
        eigenket_derivative = eigenvectors[s + 1][i] - eigenvectors[s - 1][i]
        identity_matrix = identity(2**parameters.NUM_QUBITS)
        braket_derivative = identity_matrix.matrix_element(eigenbra_derivative, eigenket_derivative)/(4*parameters.TIME_STEP)
        return braket_derivative
    
    def diabatic_tensor_components(s,i,j,k,l):
        if s == 0 or s == len(parameters.ANNEALING_PARAMETER)-1:
            return 0
        else:
            return helper.delta(i,k)*basis_dragging_matrix_components(s,l,j) + helper.delta(j,l)*basis_dragging_matrix_components(s,i,k)
    
    return [[[[[diabatic_tensor_components(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in parameters.ANNEALING_PARAMETER]


def frequency_tensor(frequency_matrix):
    states = range(len(frequency_matrix[0]))
    def frequency_tensor_components(s,i,j,k,l):
        return helper.delta(j,l)*helper.delta(i,k)*frequency_matrix[s][i][j]
    return [[[[[frequency_tensor_components(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in parameters.ANNEALING_PARAMETER]



def linbladian:

	# Initialize and diagonalize D Wave Hamiltonians
    hamiltonian = [primitive.d_wave_hamiltonian(s) for s in parameters.ANNEALING_PARAMETER]
    eigenvalues = [hamiltonian[s].eigenenergies() for s in parameters.ANNEALING_PARAMETER]
    eigenvectors = [hamiltonian[s].eigenstates() for s in parameters.ANNEALING_PARAMETER]
    idx = [eigenvalues[s].argsort()[::-1] for s in parameters.ANNEALING_PARAMETER]
    sorted_eigenvalues = [eigenvalues[s][idx[s]] for s in parameters.ANNEALING_PARAMETER]
    sorted_eigenvectors = [eigenvectors[s][:,idx[s]] for s in parameters.ANNEALING_PARAMETER]
    truncated_eigenvalues = eigenvalues_sorted[::STATE_CUT_OFF]
    
    # Raw spectral data for computation of Redfield and diabatic tensors
    truncated_eigenvectors = sorted_eigenvectors[::STATE_CUT_OFF]
    frequency_matrix = [[sorted_eigenvalues[i]-sorted_eigenvalues[j] for i in range(2**parameters.NUM_QUBITS)] for j in range(2**parameters.NUM_QUBITS)]
    truncated_frequency_matrix = frequency_matrix[::parameters.NUM_STATES_CUTOFF][::parameters.NUM_STATES_CUTOFF]
    
    # non-interacting term
    noninteracting_term = -1j*frequency_tensor(truncated_frequency_matrix)

    # Redfield term
    redfield_term = redfield.redfield_tensor(sorted_eigenvectors, frequency_matrix)
 
    # diabatic term
    diabatic_term = diabatic_tensor(truncated_eigenvectors)
    
    return noninteracting_term - (redfield_term - diabatic_term)
    



# 
