"""
Original Developer: David Roberts
Purpose of Module: To create computational primitive (F-model) from input parameters
Last Modified: 5/30/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
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


#Define D Wave Hamiltonian
def ising_coupling(qubit_index):
    if qubit_index >= parameters.NUM_QUBITS - 1:
        print "Error. coupler is over-indexed"
    else:
        if parameters.NUM_QUBITS % 2 == 0:
            if qubit_index in [parameters.NUM_QUBITS/2 - 1, parameters.NUM_QUBITS/2]:
                return -parameters.J
            else:
                return -parameters.K
        else:
            if qubit_index in [parameters.NUM_QUBITS/2 - 1, parameters.NUM_QUBITS/2]:
                return -parameters.J
            else:
                return -parameters.K
        

# Defines the computational primitive; a 1D spin glass
def f_model_problem_hamiltonian:
	qubits = range(parameters.NUM_QUBITS)
	bulk_qubits = range(parameters.NUM_QUBITS-1)

	def bulk_ising_term(qubit_index): 
		return ising_coupling(qubit_index)*helper.Z(qubit_index)*helper.Z(qubit_index + 1)

    boundary_ising_term = parameters.I*helper.Z(qubits[-1])*helper.Z(qubits[0])
    return sum([bulk_ising_term(qubit_index) for qubit_index in bulk_qubits]) + boundary_term


def driver_hamiltonian:
	qubits = range(parameters.NUM_QUBITS)
	return sum([helper.X(qubit_index) for qubit_index in qubits])


def d_wave_hamiltonian(s):
    qubits = range(parameters.NUM_QUBITS)

    def A(s):
    	return parameters.DRIVER_COEFFICIENT[s]
    def B(s):
    	return parameters.PROBLEM_COEFFICIENT[s]

    rescaled_driver_hamiltonian = (A(s)/2)*driver_hamiltonian
    rescaled_f_model_problem_hamiltonian = (B(s)/2)*f_model_problem_hamiltonian

    return 2*np.pi*(rescaled_f_model_problem_hamiltonian + rescaled_driver_hamiltonian)
                    
