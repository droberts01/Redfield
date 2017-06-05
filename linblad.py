"""
Original Developer: David Roberts
Purpose of Module: exports linblad.csv according to Redfield approximation
Last Modified: 6/5/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
from qutip import *
import numpy

# Custom Modules:
import parameters
import redfield

# Load parameters
list_of_s = parameters.S_VALUES
list_of_t = parameters.LIST_OF_TIMES

# Program runs here:
list_of_hamiltonians = map(primitive.hamiltonian, list_of_s)
list_of_eigenvalues = map(np.linalg.eigvals, list_of_hamiltonians)
list_of_eigenvalues = map(helper.truncate, list_of_eigenvalues)
list_of_frequency_tensors = map(redfield.compute_frequency_tensor, list_of_eigenvalues)
list_of_redfield_output = map(redfield.compute_redfield_tensor, list(zip(list_of_s, list_of_hamiltonians)))
list_of_redfield_tensors, list_of_eigenstates = zip(list_of_redfield_output)
list_of_eigenstates = map(helper.truncate, list_of_eigenstates)
list_of_diabatic_tensors = redfield.compute_diabatic_tensors(list_of_eigenstates)
list_of_linblads =  map(get_sum_tensors, zip(list_of_frequency_tensors, list_of_redfield_tensors, list_of_diabatic_tensors)))
list_of_linblads_csv = np.array([linblad.array for for linblad in list_of_linblads])

# Export array as csv
np.savetxt("linblad.csv", list_of_linblads_csv, delimiter=",")



