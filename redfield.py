"""
Original Developer: David Roberts
Purpose of Module: provides functions to construct the Linbladian in linblad.py,
according to Redfield approximation
Last Modified: 6/5/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
from qutip import *
import numpy as np


# Custom Modules:
import parameters
import helper



# 3. Define the Redfield tensor
# The Redfield tensor R_ijkl is a term in the Linbladian which encodes the effects of the environment
# on the system treated within the Redfield formalism.

num_qubits = parameters.NUM_QUBITS
qubits = range(num_qubits)


def compute_redfield_tensor(s, hamiltonian):
	def redfield_spectral_density(qubit):
		return spectral_density_function(s)
	redfield_tensor_Qobj, eigenstates = bloch_redfield_tensor(hamiltonian, map(helper.Z, qubits), 
																map(redfield_spectral_density, qubits))
	redfield_tensor_reals = np.real(redfield_tensor_Qobj.full())
	def compact_tensor_components(multi_index):
		I, J = multi_index
		return redfield_tensor_reals[I][J]
	redfield_compact_tensor = helper.Compact_Tensor(compact_tensor_components)
	redfield_tensor = get_tensor_from_compact_tensor(redfield_compact_tensor)
	return redfield_tensor


system_temperature = parameters.OPERATING_TEMPERATURE
bath_cutoff_time = parameters.BATH_CUTOFF_TIME
boltzmann_constant = parameters.KB
hbar = parameters.HBAR
bath_coupling = parameters.BATH_COUPLING


def spectral_density_function(s):
	def spectral_density(frequency):
		numerator = hbar**2 * parameters.bath_coupling(s) * frequency * np.exp(-abs(frequency)*bath_cutoff_time)
		demoninator = 1 - np.exp(-hbar * frequency / (boltzmann_constant * system_temperature))
		return numerator/denominator
	return spectral_density



# 2. Define the frequency tensor
# The frequency tensor defines the non-interacting dynamics of the system.

def compute_frequency_tensor(eigenvalues):
	frequency_matrix = [[Ei - Ej for Ei in eigenvalues] for Ej in eigenvalues]
	def tensor_components(multi_index):
		i, j, k, l = multi_index
		return -1j * frequency_matrix[i][j] * helper.delta(i, k) * helper.delta(j, l)
	return helper.Tensor(tensor_components)



# 3. Define the diabatic tensor
# The diabatic tensor M_ijkl is a term in the Linbladian which encodes diabatic transitions in the system.

num_eigenstates = parameters.NUM_STATES
eigenstates = range(num_eigenstates)
list_of_t = parameters.LIST_OF_TIMES
list_of_time_indices = range(len(list_of_t))


def compute_diabatic_tensors(eigenstates):
	def braket(i, j, time_index):
		bra = eigenstates[time_index + 1][i] + eigenstates[time_index - 1][i]
		ket = eigenstates[time_index + 1][j] - eigenstates[time_index - 1][j]
		time_step = list_of_t[time_index + 1] - list_of_t[time_index]
		identity_matrix = identity(num_eigenstates)
		return identity_matrix.matrix_element(bra, ket)/(4 * time_step)

	def tensor_component_function(time_index):
		if time_index == 0:
			def tensor_components(multi_index):
				i, j, k, l = multi_index
				return 0
		else:
			def tensor_components(multi_index):
				i, j, k, l = multi_index
				return helper.delta(i, k)  *braket(l, j, time_index) + helper.delta(j, l) * braket(i, k, time_index)
		return tensor_components

	list_of_diabatic_tensor_component_functions = map(tensor_component_function, list_of_time_indices)
	return map(helper.Tensor, list_of_diabatic_tensor_component_functions)





