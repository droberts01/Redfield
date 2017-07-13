#---------------
# redfield.py 
# parallel-computes linblad for the Redfield quantum master equation

import json
from pathlib import Path
import numpy as np
from scipy import linalg
import multiprocessing
import tqdm
from time import time
from functools import partial
# from matplotlib import pyplot

import meta
import meta_functions


def generate_linblad(args):
	index, svals, tvals, bad_svals, Nc, N, decoherence, H_args = args
	s_now = svals[index]
	t_now = tvals[index]
	hamiltonian_now = meta_functions.generate_hamiltonian(s_now, H_args)
	evals_now, ekets_now = linalg.eigh(hamiltonian_now)

	if not decoherence:
		if index in bad_svals:
			linblad = O(evals_now[:Nc], Nc)
		elif index not in bad_svals:
			s_before, s_after = [svals[index - 1], svals[index + 1]]
			t_before, t_after = [tvals[index - 1], tvals[index + 1]]
			dt = t_after - t_now 
			hamiltonian_before = meta_functions.generate_hamiltonian(
															s_before, H_args)
			hamiltonian_after = meta_functions.generate_hamiltonian(
															s_after, H_args)
			linblad = O(evals_now[:Nc], Nc) +\
						M(hamiltonian_before, hamiltonian_after, Nc, dt)
	elif decoherence:
		if index in bad_svals:
			linblad = O(evals_now[:Nc], Nc) - R(evals_now,
												 ekets_now, s_now, N, Nc)
		elif index not in bad_svals: 
			s_before, s_after = [svals[index - 1], svals[index + 1]]
			t_before, t_after = [tvals[index - 1], tvals[index + 1]]
			dt = t_after - t_now 
			hamiltonian_before = meta_functions.generate_hamiltonian(
														s_before, H_args)
			hamiltonian_after = meta_functions.generate_hamiltonian(
														s_after, H_args)
			linblad = O(evals_now[:Nc], Nc) +\
						M(hamiltonian_before, hamiltonian_after, Nc, dt) -\
						R(evals_now, ekets_now, s_now, N, Nc)

	return linblad


def O(evals, Nc):
	def components(i, j, k, l):
		return meta_functions.delta(i, k) * meta_functions.delta(j, l) *\
				(evals[i] - evals[j])
	return np.array([[[[- 1j * components(i, j, k, l) for l in range(Nc)] 
												for k in range(Nc)]
												for j in range(Nc)]
												for i in range(Nc)])

def R(evals, ekets, sval, N, Nc):
	unitary = ekets
	W = evals[:,np.newaxis] - evals[np.newaxis,:]

	# Compute interaction operators in instantaneous eigenbasis
	start = time()
	Z = np.array([np.matmul(
					np.transpose(unitary), 
					np.matmul(meta_functions.Z(q, N), unitary)
							) for q in range(N)])
	# Z = np.array([[[np.dot(ekets[:,j], 
	# 						np.matmul(meta_functions.Z(q, N), 
	# 													ekets[:,i]
	# 									)
	# 					) for i in range(Nc)] for j in range(Nc)] for q in range(N)]
	# 			)
	end = time()
	print("computed matrix elements of Z in {} seconds.".format(end-start))

	# print(unitary[:,0])
	# print(unitary[:,1])
	# print(Z[3,:4,:4])
	# print(Z[3,:4,:4])

	# Evaluate time-dependent noise spectral density
	S = meta_functions.S(sval)

	# Compute terms Gamma^+ and Gamma^- in the Redfield tensor
	def gamma_plus(i,j,k,l):
		# print(S(W[k,i]))
		return sum(Z[:,i,k] * Z[:,j,l]) * S(W[k,i]) / 2.

	def gamma_minus(i,j,k,l):
		return sum(Z[:,i,k] * Z[:,j,l]) * S(W[l,j]) / 2.


	def components(i,j,k,l):
		output =  - gamma_plus(i,j,k,l) - gamma_minus(i,j,k,l)
		if j == l:
			output += sum([gamma_plus(n,n,k,i) for n in range(Nc)])

		if i == k:
			output += sum([gamma_minus(n,n,j,l) for n in range(Nc)])
		return output

	return np.array([[[[components(i,j,k,l) for l in range(Nc)] 
												for k in range(Nc)]
												for j in range(Nc)]
												for i in range(Nc)])


# Test of R().
# print("running...")
# s = 0.524
# N = 8
# I = 0.3
# J = 0.9
# K = 1.4
# Nc = 3
# start = time()
# hamiltonian_now = meta_functions.generate_hamiltonian(s, [I, J, K, N])
# evals_now, ekets_now = linalg.eigh(hamiltonian_now)
# end = time()
# print("generated and diagonalized hamiltonian in {} seconds.".format(end-start))

# start = time()
# Routput = R(evals_now, ekets_now, s, N, Nc)
# end = time()
# print("R() ran in {} seconds.".format(end-start))

# for i in range(Nc):
# 	for j in range(Nc):
# 		for k in range(Nc):
# 			for l in range(Nc):
# 				print([i,j,k,l])
# 				print(Routput[i,j,k,l])


def M(hamiltonian_before, hamiltonian_after, Nc, dt):
	evals_before, ekets_before = linalg.eigh(hamiltonian_before, 
												eigvals = (0, Nc - 1))
	evals_after, ekets_after = linalg.eigh(hamiltonian_after,
												eigvals = (0, Nc - 1))
	# print(sum([(ekets_after[i,0] + ekets_before[i,0])*\
	# 			(ekets_after[i,1] - ekets_before[i,1]) for i in range(len(ekets_after[:,0]))]))
	def bkt(i,j):
		bra = ekets_after[:,i] + ekets_before[:,i]
		ket = ekets_after[:,j] - ekets_before[:,j]
		return sum(bra[:] * ket[:])/ (4 * dt)

	# print(bkt(0,1))
	def components(i,j,k,l):
		return meta_functions.delta(i, k) * bkt(l, j) + meta_functions.delta(j, l) * bkt(k, i)

	return np.array([[[[components(i,j,k,l) for l in range(Nc)] 
												for k in range(Nc)]
												for j in range(Nc)]
												for i in range(Nc)])

# Test of M().
# s = 0.235
# N = 6
# I = 0.2
# J = 0.3
# K = 1.0
# Nc = 5
# H_before = meta_functions.generate_hamiltonian(s, [I, J, K, N])
# H_after = meta_functions.generate_hamiltonian(s + 0.001, [I, J, K, N])
# dt = 5 * 10**(-15)

# Moutput = M(H_before, H_after, Nc, dt)
# for i in range(Nc):
# 	for j in range(Nc):
# 		for k in range(Nc):
# 			for l in range(Nc):
# 				print([i,j,k,l])
# 				print(Moutput[i,j,k,l])


def generate_eq_dist(args):
	sval, Nc, H_args = args
	hamiltonian = meta_functions.generate_hamiltonian(sval, H_args)
	evals, ekets = linalg.eigh(hamiltonian, eigvals = (0, Nc - 1))
	evals = [e - evals[0] for e in evals]
	boltzmann_factors = [- meta.BETA * e for e in evals]
	partition_function = sum([np.exp(boltzmann_factor) 
								for boltzmann_factor in boltzmann_factors])
	eq_dist = np.zeros((Nc, Nc))
	for j in range(Nc):
		eq_dist[j,j] = np.exp(boltzmann_factors[j]) / partition_function

	return eq_dist

# tests generate_eq_
# s = 0.235
# N = 8
# I = 0.2
# J = 0.3
# K = 1.0
# Nc = 5
# H_args = [I, J, K, N]
# output = generate_eq_dist([s, Nc, H_args])

# print(output)


def generate_json(args):
	tQA, I, J, K, N, Nc, step, window_size, num_samples, decoherence = args
	decoherence = int(decoherence)
	# args for constructing F-model Hamiltonian
	H_args = [I, J, K, int(N)]

	# args pertaining to numerical implementation of the QME
	N_args = [int(Nc), [step, window_size, int(num_samples)]]
	# print(int(num_samples))
	svals, bad_svals, tvals = meta_functions.generate_discretization(
				tQA, H_args, N_args[1])


	# provide name for generated JSON and check if JSON already exists.
	prefix = '/data/tQA, I, J, K, N, Nc, step, window_size, num_samples, decoherence = '+ str(args)
	simulation_filename = meta.ROOT + prefix + '.json'
	simulation = { 
		'svals': svals,
		'linblad_real_part': 0,
		'linblad_imaginary_part': 0,
		'rho_real_part': 0,		
		'rho_imaginary_part': 0,
		'eq_dist': 0
	}


	# Program runs here:
	print("Testing to see if JSON  already exists...")
	if Path(simulation_filename).is_file():
		print ("The following file already exists. To force the calculation, comment-out the conditional in redfield.py:")
		print(simulation_filename)
	else:
		print ("Producing the following JSON...")
		print(simulation_filename)
		print("Initializing the quantum master equation with {} time steps...".format(len(svals)
			))

		# generate linblads in parallel
		pool = multiprocessing.Pool(processes = meta.CPU_COUNT)
		iterable = np.arange(len(svals))
		linblads = pool.map(generate_linblad, 
									tqdm.tqdm(
										zip(
											iterable, 
											[svals]*len(iterable),
											[tvals]*len(iterable), 
											[bad_svals]*len(iterable), 
											[int(Nc)]*len(iterable),
											[int(N)]*len(iterable),
											[int(decoherence)]*len(iterable),
											[H_args]*len(iterable)
											)
										)
									)

		# unpack components of the linblads into the JSON dict
		# print("unpacking components of the linblads into the JSON dict...")
		# start = time()
		simulation['linblad_real_part'] = np.array(
											[np.real(l) for l in linblads]).tolist()	
		simulation['linblad_imaginary_part'] = np.array(
											[np.imag(l) for l in linblads]).tolist()	
		# end = time()
		# print("finished downloading linblads into JSON in {} seconds.".format(end-start))

		print("generating reference equilibrium distributions...")
		eq_dist = pool.map(generate_eq_dist,
									tqdm.tqdm(
										zip(
											svals,
											[int(Nc)]*len(svals),
											[H_args]*len(svals)
											)
										)
									)
		simulation['eq_dist'] = np.array(eq_dist).tolist()
		# end = time()
		# print("finished downloading equilibium distributions into JSON in {} seconds.".format(end-start))
	
		# store the JSON dict in a JSON file
		print("uploading data into the JSON file...")
		start = time()
		with open(simulation_filename, 'w') as file:
			json.dump(simulation, file)
		end = time()
		print("Process complete. uploaded new JSON file in {} seconds.".format(end-start))


