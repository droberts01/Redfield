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


# Load global variables

# Parameters of Computational Primitive (F-model) 
I = parameters.I
J = parameters.J
K = parameters.K
NUM_QUBITS = parameters.NUM_QUBITS
ANNEALING_TIME = 5*(10**(-6)))
OPERATING_TEMPERATURE = 15.5*(10**(-3))
NUM_STATES = parameters.NUM_STATES
NUM_STATES_CUTOFF = parameters.NUM_STATES_CUTOFF
DRIVER_COEFFICIENT = parameters.DRIVER_COEFFICIENT
PROBLEM_COEFFICIENT = parameters.PROBLEM_COEFFICIENT



#Define D Wave Hamiltonian
def ising_coupling(qubit_index):
    if qubit_index >= NUM_QUBITS - 1:
        print "Error. coupler is over-indexed"
    else:
        if NUM_QUBITS % 2 == 0:
            if qubit_index in [NUM_QUBITS/2 - 1, NUM_QUBITS/2]:
                return -J
            else:
                return -K
        else:
            if qubit_index in [NUM_QUBITS/2 - 1, NUM_QUBITS/2]:
                return -J
            else:
                return -K
        

# Defines the computational primitive; a 1D spin glass
def f_model_problem_hamiltonian:
	qubits = range(NUM_QUBITS)
	bulk_qubits = range(NUM_QUBITS-1)

	def bulk_ising_term(qubit_index): 
		return ising_coupling(qubit_index)*helper.Z(qubit_index)*helper.Z(qubit_index + 1)

    boundary_ising_term = I*helper.Z(qubits[-1])*helper.Z(qubits[0])
    return sum([bulk_ising_term(qubit_index) for qubit_index in bulk_qubits]) + boundary_term




def driver_hamiltonian:
	qubits = range(NUM_QUBITS)
	return sum([helper.X(qubit_index) for qubit_index in qubits])


def d_wave_hamiltonian(s):
    qubits = range(NUM_QUBITS)

    def A(s):
    	return DRIVER_COEFFICIENT[s]
    def B(s):
    	return PROBLEM_COEFFICIENT[s]

    rescaled_driver_hamiltonian = (A(s)/2)*driver_hamiltonian
    rescaled_problem_hamiltonian = (B(s)/2)*f_model_problem_hamiltonian
 
    return 2*np.pi*(rescaled_problem_hamiltonian + rescaled_driver_hamiltonian)
                    
