"""
Original Developer: David Roberts
Purpose of Module: exports linblad.csv according to Redfield approximation
Last Modified: 6/5/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
from qutip import *
import numpy as np

# Custom Modules:
import parameters
import redfield
import primitive
import helper
import time

# Load parameters
list_of_s = parameters.S_VALUES
list_of_t = parameters.LIST_OF_TIMES


globalstart = time.time()

# Program runs here:
print ("computing hamiltonians...")
start = time.time() 
list_of_hamiltonians_Qobj = map(primitive.hamiltonian, list_of_s)
list_of_hamiltonians = [np.real(hamiltonian.full()) for hamiltonian in list_of_hamiltonians_Qobj]
end = time.time()
print ("computed hamiltonians in {} seconds.".format(end-start))


print ("computing redfield tensors...")
start = time.time() 
list_of_redfield_output = map(redfield.compute_redfield_tensor, zip(list_of_s, list_of_hamiltonians_Qobj))
print ("unpacking redfield output...")
list_of_redfield_tensors, list_of_eigenstates = zip(*list_of_redfield_output)
end = time.time()
print ("computed redfield tensors in {} seconds.".format(end-start))


print ("truncacting eigenstates...")
list_of_eigenstates = map(helper.truncate, list_of_eigenstates)
start = time.time() 
print ("computing diabatic tensors...")
list_of_diabatic_tensors = redfield.compute_diabatic_tensors(list_of_eigenstates)
end = time.time()
print ("computed diabatic tensors in {} seconds.".format(end-start))


print ("constructing Linblads...")
list_of_linblads = map(helper.get_sum_tensors, zip(list_of_redfield_tensors, list_of_diabatic_tensors))
list_of_compact_linblads = map(helper.get_compact_compact_tensor_from_tensor, list_of_linblads)
list_of_compact_linblads_data = [linblad.array for linblad in list_of_compact_linblads]
list_of_compact_linblads_csv = np.array([linblad.array for linblad in list_of_compact_linblads_data])

# Export array as csv
np.savetxt("linblad.csv", list_of_linblads_csv, delimiter=",")

globalend = time.time()
print ("linblad.py complete. process took {} seconds.".format(globalend-globalstart))




