"""
Original Developer: David Roberts
Purpose of Module: Common quantum helper functions
Last Modified: 5/30/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
import numpy

# Custom Modules:
import parameters


#Define Pauli matrices

def X(acting_qubit):
    qubits = range(parameters.NUM_QUBITS)
    if acting_qubit >= parameters.NUM_QUBITS:
        print "Error. Pauli matrix over-indexed"
    else:
        def X_tensor(acting_qubit, qubit):
            if qubit == acting_qubit:
                return sigmax()
            else: 
                return identity(2)
        return tensor([X_tensor(acting_qubit, qubit) for qubit in qubits])

def Y(acting_qubit):
    qubits = range(parameters.NUM_QUBITS)
    if acting_qubit >= parameters.NUM_QUBITS:
        print "Error. Pauli matrix over-indexed"
    else:
        def Y_tensor(acting_qubit, qubit):
            if qubit == acting_qubit:
                return sigmay()
            else: 
                return identity(2)
        return tensor([Y_tensor(acting_qubit, qubit) for qubit in qubits])


def Z(acting_qubit):
    qubits = range(parameters.NUM_QUBITS)
    if acting_qubit >= parameters.NUM_QUBITS:
        print "Error. Pauli matrix over-indexed"
    else:
        def Z_tensor(acting_qubit, qubit):
            if qubit == acting_qubit:
                return sigmaz()
            else: 
                return identity(2)
        return tensor([Z_tensor(acting_qubit, qubit) for qubit in qubits])


# Define kronecker delta
def delta(i,j):
    if i == j:
        return 1
    else:
        return 0