"""
Original Developer: David Roberts
Purpose of Module: To create computational primitive used by Google QuAIL from input parameters
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
# Parameters of Google Computational Primitive
h1 = parameters.h1
h2 = parameters.h2
J = parameters.J
NUM_QUBITS = parameters.NUM_QUBITS
ANNEALING_TIME = parameters.ANNEALING_TIME
OPERATING_TEMPERATURE = parameters.OPERATING_TEMPERATURE
NUM_STATES = parameters.NUM_STATES
NUM_STATES_CUTOFF = parameters.NUM_STATES_CUTOFF
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
    rescaled_problem_hamiltonian = (B(s)/2)*google_probe_hamiltonian
 
    return 2*np.pi*(rescaled_problem_hamiltonian + rescaled_driver_hamiltonian)
              
