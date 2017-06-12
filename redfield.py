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
import math

# Custom Modules:
import parameters
import helper
import time


# 3. Define the Redfield tensor
# The Redfield tensor R_ijkl is a term in the Linbladian which encodes the effects of the environment
# on the system treated within the Redfield formalism.

list_of_s = parameters.S_VALUES
num_qubits = parameters.NUM_QUBITS
qubits = range(num_qubits)
num_states = parameters.NUM_STATES
states = range(num_states)
num_truncated_states = parameters.NUM_STATES_CUTOFF
truncated_states = range(num_truncated_states)
num_compact_truncated_states = num_truncated_states**2
compact_truncated_states = range(num_compact_truncated_states)

def compute_redfield_tensor(args):
	s, hamiltonian = args
	s_int = max(int(round(len(list_of_s) * s)) - 1, 0)
	def redfield_spectral_density(qubit):
		return Jw(s_int)
	start = time.time()
	redfield_tensor_Qobj, eigenstates = bloch_redfield_tensor(hamiltonian, map(helper.Z, qubits), 
																map(redfield_spectral_density, qubits), use_secular=False)

	end = time.time()
	print ("computed the raw redfield tensor with QuTip in {} seconds.".format(end-start))
	start = time.time()
	redfield_tensor_reals = np.real(redfield_tensor_Qobj.full())
	eigenstates = np.array([np.real(eigenstate.full()) for eigenstate in eigenstates])
	eigenstates = np.array([eigenstate[:,0] for eigenstate in eigenstates])
	eigenstates = np.transpose(eigenstates)
	end = time.time()
	print ("unpacked the redfield tensor in {} seconds.".format(end-start))
	def compact_tensor_components(multi_index):
		I, J = multi_index
		return redfield_tensor_reals[I][J]
	print ("reinitializing the tensor...")
	redfield_compact_tensor = np.array([[compact_tensor_components([I,J]) for J in compact_truncated_states]
																			for I in compact_truncated_states])
	print ("reshaping the tensor...")
	redfield_tensor = helper.get_tensor_from_compact_tensor(redfield_compact_tensor)
	print ("finished running compute_redfield_tensor(). Starting a new process...")
	return [redfield_tensor, eigenstates]




def compute_redfield_tensor_v2(args):
	s, hamiltonian = args
	s_int = max(int(round(len(list_of_s) * s)) - 1, 0)
	start = time.time()
	eigenvalues, eigenstates = np.linalg.eigh(hamiltonian)
	end = time.time()
	print ("spectral decomposition of hamiltonian computed in {} seconds.".format(end-start))
	system_part = compute_system_part(eigenvalues)
	dissipative_part = compute_dissipative_part(eigenstates, eigenvalues, s_int)
	redfield_tensor = helper.get_sum_tensors([system_part, dissipative_part]) #make sure sign is neg. on dissipative part
	return [redfield_tensor, eigenstates]


def compute_system_part(eigenvalues):
	def tensor_components(i,j,k,l):
		return -1j * helper.delta(i,k) * helper.delta(j,l) * (eigenvalues[i] - eigenvalues[j])
	return np.array([[[[tensor_components(i,j,k,l) for l in truncated_states]
													for k in truncated_states]
														for j in truncated_states]
															for i in truncated_states])


def compute_dissipative_part(eigenstates, eigenvalues, s_int):
	U = eigenstates
	W = eigenvalues[:,np.newaxis] - eigenvalues[np.newaxis,:]
	start = time.time()
	Z = np.array([np.matmul(np.transpose(U), np.matmul(np.real(helper.Z(qubit).full()), U)) 
																			for qubit in qubits])
	end = time.time()
	print("computed Z in {} seconds.".format(end-start))

	Jws = Jw(s_int)
	def G_plus(i,j,k,l):
		return Jws(-W[j,l]) * sum(Z[:,i,k] * Z[:,j,l])/2
	def G_minus(i,j,k,l):
		return Jws(W[i,k]) * sum(Z[:,i,k] * Z[:,j,l])/2

	def tensor_components(i,j,k,l):
		sum_G_plus_innk = sum([G_plus(i,n,n,k) for n in truncated_states])
		n = num_truncated_states
		while np.abs(W[n,k]) < 10**11 and n < num_states - 1:
			sum_G_plus_innk = sum_G_plus_innk + G_plus(i,n,n,k)
			n = n + 1

		sum_G_minus_jnnl = sum([G_plus(j,n,n,l) for n in truncated_states])
		n = num_truncated_states
		while np.abs(W[n,j]) < 10**11 and n < num_states - 1:
			sum_G_minus_jnnl = sum_G_minus_jnnl + G_minus(j,n,n,l)
			n = n + 1

		part_one = helper.delta(j,l) * sum_G_plus_innk  + helper.delta(i,k) * sum_G_minus_jnnl
		part_two = - G_plus(i,j,k,l) - G_minus(i,j,k,l)
		return part_one + part_two
	return np.array([[[[-tensor_components(i,j,k,l) for l in truncated_states]
													for k in truncated_states]
														for j in truncated_states]
															for i in truncated_states])




system_temperature = parameters.OPERATING_TEMPERATURE
bath_cutoff_time = parameters.BATH_CUTOFF_TIME
boltzmann_constant = parameters.KB
hbar = parameters.HBAR
bath_coupling = parameters.BATH_COUPLING


spectral_density_coefficient = [hbar**2 * bath_coupling(s) for s in list_of_s]
exponential_coefficient = hbar / (boltzmann_constant * system_temperature)

def Jw(s_int):
	def spectral_density(frequency):
		numerator = spectral_density_coefficient[s_int] * frequency 
		# Analytical continuation of spectral density to the real axis:
		if abs(frequency) > 10**11: #Prevent numerical overflow of np.exp() - 100 GHz cutoff
			return 0.0
		elif frequency != 0:
			denominator = 1 - np.exp(-exponential_coefficient*frequency)
			return numerator / denominator
		else:
			return spectral_density_coefficient[s_int] / exponential_coefficient
	return spectral_density


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
		return sum(bra[:] * ket[:])/(4 * time_step)

	def tensor_component_function(time_index):
		if time_index == 0 or time_index >= len(list_of_t)- 1:
			def tensor_components(multi_index):
				i, j, k, l = multi_index
				return 0
		else:
			def tensor_components(multi_index):
				i, j, k, l = multi_index
				return helper.delta(i, k)  *braket(l, j, time_index) + helper.delta(j, l) * braket(i, k, time_index)
		return tensor_components

	list_of_diabatic_component_functions = map(tensor_component_function, list_of_time_indices)
	return [np.array([[[[list_of_diabatic_component_functions[time_index]([i,j,k,l]) 
																			for l in truncated_states]
																				for k in truncated_states]
																					for j in truncated_states]
																						for i in truncated_states])
																							for time_index in list_of_time_indices]


