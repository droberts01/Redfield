"""
Original Developer: David Roberts
Purpose of Module: Common quantum helper functions
Last Modified: 6/5/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
from qutip import *
import numpy as np

# Custom Modules:
import parameters


# Load parameters
num_qubits = parameters.NUM_QUBITS
qubits = range(num_qubits)




# Rank-four tensor.

num_eigenstates = parameters.NUM_STATES_CUTOFF
eigenstates = range(num_eigenstates)

class Tensor:
    def __init__(self, components):
        self.array = np.array([[[[components([i,j,k,l]) for l in eigenstates] 
                                                                for k in eigenstates] 
                                                                    for j in eigenstates] 
                                                                        for i in eigenstates])

        self.components = components




# Vectorized tensor, a.k.a. superoperator


num_compact_eigenstates = num_eigenstates**2
compact_eigenstates = range(num_compact_eigenstates)


class Compact_Tensor:
    def __init__(self, components):
        self.array = np.array([[components([I,J]) for J in compact_eigenstates] 
                                                                for I in compact_eigenstates])
        self.components = components



# Completely vectorized tensor for .csv export.

num_compact_compact_eigenstates = num_compact_eigenstates**2
compact_compact_eigenstates = range(num_compact_compact_eigenstates)


class Compact_Compact_Tensor:
    def __init__(self, components):
        self.array = np.array([components(I) for I in compact_compact_eigenstates])
        self.components = components


def init_compact_compact_tensor_from_array(array):
    def compact_compact_components(index):
        return array[index]
    return Compact_Compact_Tensor(compact_compact_components)



def get_sum_tensors(list_of_tensors):
    output = np.zeros((num_eigenstates, num_eigenstates, num_eigenstates, num_eigenstates), dtype = float)
    for tensor in list_of_tensors:
        output = np.add(output, tensor)
    return output

def get_sum_compact_tensors(list_of_compact_tensors):
    output = np.zeros((num_compact_eigenstates, num_compact_eigenstates), dtype = float)
    for tensor in list_of_compact_tensors:
        output = np.add(output, tensor)
    return output


def get_sum_compact_compact_tensors(list_of_compact_compact_tensors):
    output = np.zeros((num_compact_compact_eigenstates), dtype = float)
    for tensor in list_of_compact_compact_tensors:
        output = np.add(output, tensor)
    return output



# Utilizes formula R_ijkl = R_IJ, where I = Ni + j,  J = Nk + l
def get_tensor_from_compact_tensor(args):
    compact_tensor, dim = args
    def components(multi_index):
        i, j, k, l = multi_index
        I = dim * i + j
        J = dim * k + l
        return compact_tensor[I, J]
    return np.array([[[[components([i,j,k,l]) for l in range(dim)]
                                                for k in range(dim)]
                                                    for j in range(dim)]
                                                        for i in range(dim)])


# Utilizes formula R_I = R_ijkl, where I = i N^3 + j N^2 + k N + l
def get_compact_compact_tensor_from_tensor(args):
    tensor, dim = args
    def compact_compact_components(index):
        i,remainder = divmod(index, dim**3)
        j,remainder = divmod(remainder, dim**2)
        k,l = divmod(remainder, dim)
        return tensor[i, j, k, l]
    return np.array([compact_compact_components(I) for I in range(dim**4)])

def get_compact_compact_tensor_from_compact_tensor(args):
    compact_tensor, dim = args
    def compact_compact_components(index):
        I, J = divmod(index, dim**2)
        return compact_tensor[I, J]
    return np.array([compact_compact_components(I) for I in range(dim**4)])


def get_tensor_from_compact_compact_tensor(args):
    compact_compact_tensor, dim = args
    def components(multi_index):
        i, j, k, l = multi_index
        return compact_compact_tensor[dim**3 * i + 
                                            dim**2 * j + dim * k + l]
    return np.array([[[[components([i,j,k,l]) for l in range(dim)]
                                                for k in range(dim)]
                                                    for j in range(dim)]
                                                        for i in range(dim)])


def get_compact_tensor_from_compact_compact_tensor(args):
    compact_compact_tensor, dim = args
    def compact_components(multi_index):
        I, J = multi_index
        return compact_compact_tensor[dim**2 * I + J]
    return np.array([[compact_components([I,J]) for J in range(dim**2)]
                                                    for I in range(dim**2)])






# Define Pauli matrices
def X(acting_qubit):
    qubits = range(num_qubits)
    if acting_qubit >= num_qubits:
        print "Error. Pauli matrix over-indexed"
    else:
        def X_tensor(acting_qubit, qubit):
            if qubit == acting_qubit:
                return sigmax()
            else: 
                return identity(2)
        return tensor([X_tensor(acting_qubit, qubit) for qubit in qubits])

def Y(acting_qubit):
    qubits = range(num_qubits)
    if acting_qubit >= num_qubits:
        print "Error. Pauli matrix over-indexed"
    else:
        def Y_tensor(acting_qubit, qubit):
            if qubit == acting_qubit:
                return sigmay()
            else: 
                return identity(2)
        return tensor([Y_tensor(acting_qubit, qubit) for qubit in qubits])


def Z(acting_qubit):
    qubits = range(num_qubits)
    if acting_qubit >= num_qubits:
        print "Error. Pauli matrix over-indexed"
    else:
        def Z_tensor(acting_qubit, qubit):
            if qubit == acting_qubit:
                return sigmaz()
            else: 
                return identity(2)
        return tensor([Z_tensor(acting_qubit, qubit) for qubit in qubits])



# Helper functions
def truncate(list_of_objects):
    return list_of_objects[:parameters.NUM_STATES_CUTOFF]


# Define Kronecker delta
def delta(i,j):
    if i == j:
        return 1
    else:
        return 0
