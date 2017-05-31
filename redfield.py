"""
Original Developer: David Roberts
Purpose of Module: outputs Redfield tensor
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
import helper

# Load global variables
NUM_QUBITS = parameters.NUM_QUBITS
NUM_STATES = parameters.NUM_STATES
SPECTRAL_DENSITY = parameters.SPECTRAL_DENSITY


def sigmaz_matrix_element(s,qubit_index,i,j, sorted_eigenvectors):
            helper.Z(qubit_index).matrix_element(sorted_eigenvectors[s][i],sorted_eigenvectors[s][j])


def redfield_coefficient_tensor(sorted_eigenvectors):
 
    def redfield_coefficient_tensor_components(s,i,j,k,l):
        return sum([sigmaz_matrix_element(s,qubit,i,k, sorted_eigenvectors)*sigmaz_matrix_element(s,qubit,j,l, sorted_eigenvectors) for qubit in range(NUM_QUBITS)])
 
    return [[[[[redfield_coefficient_tensor_components(s,i,j,k,l) for l in range(NUM_STATES)] for k in range(NUM_STATES)] for j in range(NUM_STATES)] for i in range(NUM_STATES)] for s in ANNEALING_PARAMETER]


def redfield_tensor(sorted_eigenvectors, frequency_matrix):
    states = range(NUM_STATES)
    qubits = range(NUM_QUBITS)

    coefficient_tensor = redfield_coefficient_tensor(sorted_eigenvectors)
    
    def gamma_plus(s,i,j,k,l):
        return SPECTRAL_DENSITY(-frequency_matrix[s][k][l])*coefficient_tensor[s][i][j][k][l]
    def gamma_minus(s,i,j,k,l):
        return SPECTRAL_DENSITY(frequency_matrix[s][i][j])*coefficient_tensor[s][i][j][k][l]

    def redfield_tensor_term_no1_components(s,i,j,k,l):
        return sum([helper.delta(j,l)*gamma_plus(s,i,k,n,n) for n in range(2**num_qubits)])

    def redfield_tensor_term_no2_components(s,i,j,k,l):
        return sum([helper.delta(i,k)*gamma_minus(s,j,l,n,n) for n in range(2**num_qubits)])

    def redfield_tensor_term_no3_components(s,i,j,k,l):
        return gamma_plus(s,i,j,k,l)

    def redfield_tensor_term_no4_components(s,i,j,k,l):
        return gamma_minus(s,i,j,k,l)
    
    redfield_tensor_term_no1 = [[[[[redfield_tensor_term_no1_components(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]
    redfield_tensor_term_no2 = [[[[[redfield_tensor_term_no2_components(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]
    redfield_tensor_term_no3 = [[[[[redfield_tensor_term_no3_components(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]
    redfield_tensor_term_no4 = [[[[[redfield_tensor_term_no4_components(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]
    
    return redfield_tensor_term_no1 + redfield_tensor_term_no2 - redfield_tensor_term_no3 - redfield_tensor_term_no4

