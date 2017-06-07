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




# get_sum_of_functions() (JW)
def get_sum_of_functions(list_of_functions):
     def sum_of_functions(inputs):
         outputs = [function(inputs) for function in list_of_functions]
         sum_outputs = sum(outputs)
         return sum_outputs
     return sum_of_functions



def get_sum_tensors(list_of_tensors):
    list_of_components = [tensor.components for tensor in list_of_tensors]
    sum_of_components = get_sum_of_functions(list_of_components)
    return Tensor(sum_of_components)


def get_sum_compact_tensors(list_of_compact_tensors):
    list_of_components = [tensor.components for tensor in list_of_compact_tensors]
    sum_of_components = get_sum_of_functions(list_of_components)
    return Compact_Tensor(sum_of_components)


def get_sum_compact_compact_tensors(list_of_compact_compact_tensors):
    list_of_components = [tensor.components for tensor in list_of_compact_compact_tensors]
    sum_of_components = get_sum_of_functions(list_of_components)
    return Compact_Compact_Tensor(sum_of_components)



# Utilizes formula R_ijkl = R_IJ, where I = Ni + j,  J = Nk + l
def get_tensor_from_compact_tensor(compact_tensor):
    compact_components = compact_tensor.components 
    def components(multi_index):
        i, j, k, l = multi_index
        I = num_eigenstates * i + j
        J = num_eigenstates * k + l
        return compact_components([I, J])
    return Tensor(components)


# Utilizes formula R_I = R_ijkl, where I = i N^3 + j N^2 + k N + l
def get_compact_compact_tensor_from_tensor(tensor):
    components = tensor.components
    def compact_compact_components(index):
        i,remainder = divmod(index, num_eigenstates**3)
        j,remainder = divmod(remainder, num_eigenstates**2)
        k,l = divmod(remainder, num_eigenstates)
        return components([i, j, k, l])
    return Compact_Compact_Tensor(components)


def get_tensor_from_compact_compact_tensor(compact_compact_tensor):
    compact_compact_components = compact_compact_tensor.components
    def components(multi_index):
        i, j, k, l = multi_index
        return compact_compact_components((num_eigenstates**3 * i + 
                                            num_eigenstates**2 * j + num_eigenstates * k + l))
    return Tensor(components)


def get_compact_tensor_from_compact_compact_tensor(compact_compact_tensor):
    compact_compact_components = compact_compact_tensor.components
    def compact_components(multi_index):
        I, J = multi_index
        return compact_compact_components((num_eigenstates**2 * I + J))
    return Compact_Tensor(compact_components)

def get_compact_compact_tensor_from_compact_tensor(compact_tensor):
    compact_components = compact_tensor.components
    def compact_compact_components(index):
        I, J = divmod(index, num_eigenstates**2)
        return compact_components([I, J])
    return Compact_Compact_Tensor(compact_compact_components)




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
