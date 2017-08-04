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


# generate_linblad(999, )

def generate_linblad(args):
	index, svals, tvals, bad_svals, Nc, N, decoherence, LF_noise, evals, eket_3window = args
	s_now = svals[index]
	t_now = tvals[index]
	# print(index)

	if not decoherence:
		if index in bad_svals:
			evals_now, ekets_now = [evals, eket_3window]		
			linblad = O(evals_now[:Nc], Nc)
		elif index not in bad_svals:
			# print(bad_svals)
			# print(eket_3window)
			s_before, s_after = [svals[index - 1], svals[index + 1]]
			t_before, t_after = [tvals[index - 1], tvals[index + 1]]
			dt = t_after - t_now 
			evals_now, ekets_now = [evals, eket_3window[1]]		
			ekets_before = eket_3window[0]
			ekets_after = eket_3window[2]
			linblad = O(evals_now[:Nc], Nc) +\
						M(ekets_before[:,:Nc], ekets_after[:,:Nc], Nc, dt)
	elif decoherence:
		if index in bad_svals:
			evals_now, ekets_now = [evals, eket_3window]		
			linblad = O(evals_now[:Nc], Nc) - R(evals_now, 
											ekets_now, s_now, N, Nc, LF_noise)
		elif index not in bad_svals: 
			s_before, s_after = [svals[index - 1], svals[index + 1]]
			t_before, t_after = [tvals[index - 1], tvals[index + 1]]
			dt = t_after - t_now
			evals_now, ekets_now = [evals, eket_3window[1]]		
			ekets_before = eket_3window[0]
			ekets_after = eket_3window[2]
			# print("index is {}.".format(index))
			# print ("ekets_now is:")
			# print ekets_now
			linblad = O(evals_now[:Nc], Nc) +\
						M(ekets_before[:,:Nc], ekets_after[:,:Nc], Nc, dt) -\
						R(evals_now, ekets_now, s_now, N, Nc, LF_noise)
	# print(index)
	# print(bad_svals)
	return linblad


def O(evals, Nc):
	def components(i, j, k, l):
		return meta_functions.delta(i, k) * meta_functions.delta(j, l) *\
				(evals[i] - evals[j])
	return np.array([[[[- 1j * components(i, j, k, l) for l in range(Nc)] 
												for k in range(Nc)]
												for j in range(Nc)]
												for i in range(Nc)])

def R(evals, ekets, sval, N, Nc, LF_noise):
	unitary = ekets
	# print (ekets[:,:10])
	W = evals[:,np.newaxis] - evals[np.newaxis,:]

	# Compute interaction operators in instantaneous eigenbasis
	start = time()
	Z = np.array([np.matmul(
					np.transpose(unitary), 
					np.matmul(meta_functions.Z(q, N), unitary)
							) for q in range(N)])
	end = time()
	# print("computed matrix elements of Z in {} seconds.".format(end-start))

	# Evaluate time-dependent noise spectral density
	S = meta_functions.S(sval, LF_noise)

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


# Define a smooth global O(Nc)-frame over 0 < s < 1
# def gauge_fix(ekets_before, ekets_after):
# 	for j in range(len(ekets_before)):
# 		old_diff = np.linalg.norm(ekets_before[j] - ekets_after[j])
# 		gauge_diff = np.linalg.norm(ekets_before[j] + ekets_after[j])
# 		if gauge_diff < old_diff:
# 			ekets_after[j] = - ekets_after[j]
# 	return ekets_before, ekets_after


def M(ekets_before, ekets_after, Nc, dt):
	# evals_before, ekets_before = \
	# 	meta_functions.diagonalization_func(Nc)(hamiltonian_before)
	# evals_after, ekets_after = \
	# 	meta_functions.diagonalization_func(Nc)(hamiltonian_after)

	def bkt(i,j):
		bra = ekets_after[:,i] + ekets_before[:,i]
		ket = ekets_after[:,j] - ekets_before[:,j]
		return sum(bra[:] * ket[:])/ (4 * dt)


	def components(i,j,k,l):
		return meta_functions.delta(i, k) * bkt(l, j) + meta_functions.delta(j, l) * bkt(k, i)
		# return 0

	return np.array([[[[components(i,j,k,l) for l in range(Nc)] 
												for k in range(Nc)]
												for j in range(Nc)]
												for i in range(Nc)])



def generate_eq_dist(args):
	sval, Nc, evals = args
	evals = [e - evals[0] for e in evals]
	boltzmann_factors = [- meta.BETA * e for e in evals]
	# for j in range(len(boltzmann_factors)):
	# 	print("j is {}".format(j))
	# 	print(boltzmann_factors[j])
	partition_function = sum([np.exp(boltzmann_factor) 
								for boltzmann_factor in boltzmann_factors[:Nc]])
	eq_dist = np.zeros((Nc, Nc))
	for j in range(Nc):
		eq_dist[j,j] = np.exp(boltzmann_factors[j]) / partition_function
	# for j in range(Nc):
	# 	print(j)
	# 	print(eq_dist[j,j])
	return eq_dist


def compare(args):
	ekets_before, ekets_after = args
	Nc = len(ekets_before[0,:])
	# print(ekets_before[:,0])
	def compare_individual_vector(j):
		diff = np.linalg.norm(ekets_after[:,j] - ekets_before[:,j])
		transformed_diff = np.linalg.norm(
								ekets_after[:,j] + ekets_before[:,j])
		if transformed_diff < diff:
			return -1
		else:
			return +1
	overlaps = [compare_individual_vector(j) 
								for j in range(Nc)]
	return overlaps



# holonomy starts @t=1, not t=0.
def generate_frame(evals, ekets, Nc):
	pool = multiprocessing.Pool(processes = meta.CPU_COUNT)
	list_of_overlaps = np.array(map(compare, 
		[[ekets[t], ekets[t+1]] for t in range(len(ekets)-1)]))

	holonomy = np.zeros((len(list_of_overlaps) + 1, Nc), dtype = int)
	holonomy[0] = [1]*Nc
	for t in range(len(holonomy) - 1):
		holonomy[t + 1] = [holonomy[t, j]*list_of_overlaps[t, j] 
														for j in range(Nc)]
	# for h in holonomy:
	# 	print(h)
	# x = 2 * y
	for t in range(1, len(ekets)):
		for j in range(Nc):
			ekets[t][:,j] = holonomy[t, j] * ekets[t][:,j]

	return evals, ekets



def generate_json(args):
	tQA, I, J, K, N, Nc, step, window_size, num_samples, decoherence, LF_noise, CPU = args
	decoherence = int(decoherence)
	LF_noise = int(LF_noise)
	Nc = int(Nc)
	# args for constructing F-model Hamiltonian
	H_args = [I, J, K, int(N)]

	# args pertaining to numerical implementation of the QME
	N_args = [int(Nc), [step, window_size, int(num_samples)]]
	# print(int(num_samples))
	svals, bad_svals, tvals = meta_functions.generate_discretization(
				tQA, H_args, N_args[1])


	# provide name for generated JSON and check if JSON already exists.
	prefix = 'tQA, I, J, K, N, Nc, step, window_size, num_samples, decoherence, LF_noise, CPU = '+ str(args)
	if CPU == "Darwin":
		simulation_filename = prefix + '.json'
	else:
		simulation_filename = meta.SAVE_LOCATION + prefix + '.json'
	
	simulation = { 
		'svals': svals.tolist(),
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
		# iterable = np.arange(len(svals))
		# print(svals)
		print("storing Hamiltonians in RAM...")
		hamiltonians = pool.map(meta_functions.generate_hamiltonian, 
											tqdm.tqdm(
												zip(
													svals,
													[H_args]*len(svals)
													)
												)
											)

		# H = hamiltonians[0]
		# diagonalize = partial(linalg.eigh, eigvals = (0, Nc - 1))
		# evals, ekets = diagonalize(H)
		# print(evals[0])
		# meta_functions.diagonalization_func(Nc)(hamiltonians[0])
		
		# evals, ekets = \
		# 	zip(map(meta_functions.diagonalization_func(Nc), 
		# 					hamiltonians
		# 				)
		# 		)
		print("extracting spectral data from Hamiltonians...")
		# diagonalize = meta_functions.diagonalization_func(Nc)
		data = \
			pool.map(linalg.eigh, 
							tqdm.tqdm(hamiltonians)
						)
		evals = [d[0] for d in data]
		ekets = [d[1] for d in data]
		# print(len(ekets))

		evals, ekets = generate_frame(evals, ekets, Nc)

		ekets_3windows = [ekets[0]] +\
			[ekets[j-1:j+2] for j in range(1,len(ekets)-1)] +\
												[ekets[-1]]

		for index in range(len(ekets_3windows)):
			if index in bad_svals:
				ekets_3windows[index] = ekets[index]


		# print ("ekets_3windows[280] is")
		# print (ekets_3windows[280])
		# print ("ekets_3windows[281] is")
		# print (ekets_3windows[281])
		# print ("ekets_3windows[282] is")
		# print (ekets_3windows[282])

		# print("ekets_3windows[0] is:")
		# print(ekets_3windows[0])
		# print("ekets_3windows[1] is:")
		# print(ekets_3windows[1])
		# print("we made it up to here.")


		print("generating Linblads...")
		linblads = pool.map(generate_linblad, 
									tqdm.tqdm(
										zip(
											range(len(svals)), 
											[svals]*len(svals),
											[tvals]*len(svals), 
											[bad_svals]*len(svals), 
											[int(Nc)]*len(svals),
											[int(N)]*len(svals),
											[int(decoherence)]*len(svals),
											[int(LF_noise)]*len(svals),
											evals,
											ekets_3windows
											)
										)
									)

		# unpack components of the linblads into the JSON dict
		simulation['linblad_real_part'] = np.array(
											[np.real(l) for l in linblads]).tolist()	
		simulation['linblad_imaginary_part'] = np.array(
											[np.imag(l) for l in linblads]).tolist()	



		# for k in range(int(Nc)):
		# 	for l in range(int(Nc)):
		# 		print(sum([linblads[200][j,j,k,l] for j in range(int(Nc))]))


		print("generating reference equilibrium distributions...")
		eq_dist = pool.map(generate_eq_dist,
									tqdm.tqdm(
										zip(
											svals,
											[int(Nc)]*len(svals),
											evals
											)
										)
									)
		simulation['eq_dist'] = np.array(eq_dist).tolist()

	

		# store the JSON dict in a JSON file
		print("uploading data into the JSON file...")

		start = time()
		with open(simulation_filename, 'w') as file:
			json.dump(simulation, file)
		end = time()
		print("Process complete. uploaded new JSON file in {} seconds.".format(end-start))


# Test of M().
# s = 0.235
# N = 6
# I = 0.2
# J = 0.3
# K = 1.0
# Nc = 5
# H_before = meta_functions.generate_hamiltonian(s, [I, J, K, N])
# H_after = meta_functions.generate_hamiltonian(s + 0.001, [I, J, K, N])
# dt = 5 * 10**(-9)

# Moutput = M(H_before, H_after, Nc, dt)
# # for i in range(Nc):
# # 	for j in range(Nc):
# # 		for k in range(Nc):
# # 			for l in range(Nc):
# # 				print([i,j,k,l])
# # 				print(Moutput[i,j,k,l])
# for k in range(Nc):
# 	for l in range(Nc):
# 		print([k,l])
# 		print(sum([Moutput[j,j,k,l] for j in range(Nc)]))


# # Test of R().
# print("running...")
# s = 0.235
# N = 6
# I = 0.2
# J = 0.3
# K = 1.0
# Nc = 5
# start = time()
# hamiltonian_now = meta_functions.generate_hamiltonian(s, [I, J, K, N])
# evals_now, ekets_now = linalg.eigh(hamiltonian_now)
# end = time()
# print("generated and diagonalized hamiltonian in {} seconds.".format(end-start))

# start = time()
# Routput = R(evals_now, ekets_now, s, N, Nc)
# end = time()
# print("R() ran in {} seconds.".format(end-start))

# # for i in range(Nc):
# # 	for j in range(Nc):
# # 		for k in range(Nc):
# # 			for l in range(Nc):
# # 				print([i,j,k,l])
# # 				print(Routput[i,j,k,l])
# for k in range(Nc):
# 	for l in range(Nc):
# 		print([k,l])
# 		print(sum([Routput[j,j,k,l] for j in range(Nc)]))


# tests generate_eq_dist()
# s = 0.235
# N = 8
# I = 0.2
# J = 0.3
# K = 1.0
# Nc = 5
# H_args = [I, J, K, N]
# output = generate_eq_dist([s, Nc, H_args])

# print(output)
