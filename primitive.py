"""
Original Developer: David Roberts
Purpose of Module: To create NASA and Google QuAIL computational primitives
from input parameters
Last Modified: 6/5/17
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
MODE = parameters.MODE


# Parameters of NASA Computational Primitive (F-model) 
if MODE == "NASA":
    I = parameters.I
    J = parameters.J
    K = parameters.K
    NUM_QUBITS = parameters.NUM_QUBITS
    NUM_STATES = parameters.NUM_STATES
    A = parameters.A
    B = parameters.B
else:
# Load global variables
# Parameters of Google Computational Primitive
    h1 = parameters.h1
    h2 = parameters.h2
    J = parameters.J
    NUM_QUBITS = parameters.NUM_QUBITS
    NUM_STATES = parameters.NUM_STATES
    PAIR_SET = parameters.PAIR_SET


def google_probe_hamiltonian:
    qubits = range(NUM_QUBITS)

    def onsite_hamiltonian:
        def external_field(qubit_index):
            if qubit_index < 8:
                return h1
            else:
                return h2
        def onsite_term(qubit_index):
            return external_field(qubit_index)*helper.Z(qubit_index)
        return sum([onsite_term(qubit_index) for qubit_index in qubits])


    def ising_hamiltonian:
        output = ([0]*NUM_STATES)*NUM_STATES
 
        for pair in PAIR_SET:
            output = output - J*helper.Z(pair[0])*helper.Z(pair[1])
        return output

    return onsite_hamiltonian + ising_hamiltonian



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




def hamiltonian(s):
    qubits = range(NUM_QUBITS)

    rescaled_driver_hamiltonian = (A(s)/2)*driver_hamiltonian
    if MODE == "GOOGLE":
        rescaled_problem_hamiltonian = (B(s)/2)*google_probe_hamiltonian
    else:
        rescaled_problem_hamiltonian = (B(s)/2)*f_model_problem_hamiltonian
 
    return 2*np.pi*(rescaled_problem_hamiltonian + rescaled_driver_hamiltonian)
                    
