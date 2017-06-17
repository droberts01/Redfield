
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
	redfield_compact_tensor = redfield_tensor_Qobj.full()

	dim = num_states
	redfield_tensor = helper.get_tensor_from_compact_tensor([redfield_compact_tensor, dim])
	eigenstates = np.array([eigenstate.full() for eigenstate in eigenstates])
	eigenstates = np.array([eigenstate[:,0] for eigenstate in eigenstates])
	eigenstates = np.transpose(eigenstates)

	# Truncation of R_I^J occurs here.
	truncated_redfield_tensor = redfield_tensor[:num_truncated_states,:num_truncated_states,
												:num_truncated_states,:num_truncated_states]


	# NOTE: REDFIELD TENSOR IN QUTIP SATISFIES (Rqutip)_ij^kl = (Ractual)_ji^lk.
	transposed_truncated_redfield_tensor = np.array([[[[truncated_redfield_tensor[j,i,l,k] 
																	for l in truncated_states]
																		for k in truncated_states]
																			for j in truncated_states]
																				for i in truncated_states])

	print ("finished running compute_redfield_tensor(). Starting a new process...")
	return [transposed_truncated_redfield_tensor, eigenstates]




def compute_redfield_tensor_v2(args):
	s, hamiltonian = args
	s_int = max(int(round(len(list_of_s) * s)) - 1, 0)
	start = time.time()
	eigenvalues, eigenstates = np.linalg.eigh(hamiltonian)
	end = time.time()
	# print ("spectral decomposition of hamiltonian computed in {} seconds.".format(end-start))
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
	# print("computed matrix elements of interaction operators in {} seconds.".format(end-start))

	Jws = Jw(s_int)
	def A_plus(i,j,k,l):
		return 0.5 * sum(Z[:,i,k] * Z[:,j,l]) * Jws(W[k,i])
	def A_minus(i,j,k,l):
		return 0.5 * sum(Z[:,i,k] * Z[:,j,l]) * Jws(W[l,j])

	def tensor_components(i,j,k,l):
		output =  - A_plus(i,j,k,l) - A_minus(i,j,k,l)
		if j == l:
			sum_A_plus_nnki = sum([A_plus(n,n,k,i) for n in truncated_states])
			if num_truncated_states < num_states:
				n = num_truncated_states
				while np.abs(W[k,n]) < 10**11 and n < num_states - 1:
					sum_A_plus_nnki = sum_A_plus_nnki + A_plus(n,n,k,i)
					n = n + 1

			output += sum_A_plus_nnki 
		if i == k:
			sum_A_minus_nnjl = sum([A_minus(n,n,j,l) for n in truncated_states])
			if num_truncated_states < num_states:
				n = num_truncated_states
				while np.abs(W[l,n]) < 10**11 and n < num_states - 1:
					sum_A_minus_nnjl = sum_A_minus_nnjl + A_minus(n,n,j,l)
					n = n + 1

			output += sum_A_minus_nnjl

		return output

	# start = time.time()
	dissipative_part = np.array([[[[-tensor_components(i,j,k,l) for l in truncated_states]
														for k in truncated_states]
															for j in truncated_states]
																for i in truncated_states])
	# end = time.time()
	# print("while loops ran in {} seconds. Starting a new process....".format(end-start))
	return dissipative_part


system_temperature = parameters.OPERATING_TEMPERATURE
bath_cutoff_time = parameters.BATH_CUTOFF_TIME
boltzmann_constant = parameters.KB
hbar = parameters.HBAR
bath_coupling = parameters.BATH_COUPLING


# REMOVED FACTOR OF HBAR^2
# spectral_density_coefficient = [hbar**2 * bath_coupling(s) for s in list_of_s]
spectral_density_coefficient = [ bath_coupling(s) for s in list_of_s]
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
				return helper.delta(i, k)  *braket(l, j, time_index) + helper.delta(j, l) * braket(k, i, time_index)
		return tensor_components

	return [np.array([[[[tensor_component_function(time_index)([i,j,k,l]) 
																			for l in truncated_states]
																				for k in truncated_states]
																					for j in truncated_states]
																						for i in truncated_states])
																							for time_index in list_of_time_indices]


