import numpy as np
import meta_functions
import meta
import multiprocessing
import tqdm
from time import time
import json
from pathlib import Path



def update(rho_now, L_now, L_after, dt):
	Nc2 = len(rho_now)
	I = np.identity(Nc2, dtype = complex)
	# print(I.shape)
	# print(L_after.shape)

	A_after = np.linalg.inv((I - (dt/2.)*L_after))
	A_now = (I + (dt/2.) * L_now)

	return np.matmul(A_after, np.matmul(A_now, rho_now))



def master_eq_solve(initial_condition, L, tvals, Nc):
	N_steps = len(L)
	# Initialize time-evolution as constant array
	rho = np.array([initial_condition]*(N_steps), 
										dtype=complex)


	for t in range(N_steps - 1):
		rho[t + 1] = update(rho[t], L[t], L[t + 1], 
											tvals[t + 1] - tvals[t])
		
		# Intermittently check that evolution is trace-preserving
		# if t%round((len(tvals)/20)) == 0:
		# 	print("time_index is {}".format(t))
		# 	print ("time_step is {}".format(tvals[t + 1] - tvals[t]))
		# 	print(np.transpose(rho[t - 1, np.newaxis]))
		# 	print("tr(rho[time_index]) is {}".format(sum([rho[t-1, (Nc + 1)*j]  for j in range(Nc)])))


	# Repackage density matrix back into a matrix
	rho = map(meta_functions.matrix,
									zip(
										rho, 
										[Nc]*len(rho)
										)
									)
	return rho



def solve_json(args):
	tQA, I, J, K, N, Nc, step, window_size, num_samples, decoherence, LF_noise, CPU = args
	decoherence = int(decoherence)
	LF_noise = int(LF_noise)
	# args for constructing F-model Hamiltonian
	H_args = [I, J, K, int(N)]

	# args pertaining to numerical implementation of the QME
	N_args = [int(Nc), [step, window_size, int(num_samples)]]
	svals, bad_svals, tvals = meta_functions.generate_discretization(
				tQA, H_args, N_args[1])

	# Check if JSON exists.
	prefix = 'tQA, I, J, K, N, Nc, step, window_size, num_samples, decoherence, LF_noise, CPU = '+ str(args)
	
	if CPU == "Darwin":
		simulation_filename = prefix + '.json'
	else:
		simulation_filename = meta.SAVE_LOCATION + prefix + '.json'
	
	# Program runs here:
	print("Testing to see if JSON exists...")
	if not Path(simulation_filename).is_file():
		print ("The following file does not exist. To initialize the file, run generate.py:")
		print(simulation_filename)
	else:
		print ("Solving the quantum master equation in the following JSON...")
		print(simulation_filename)
		print("loading JSON...")
		simulation = json.loads(open(simulation_filename).read())
		print("Finished loading JSON. Initializing the quantum master equation with {} time steps...".format(len(svals)
			))

		initial_condition = np.array([0.]*(int(Nc**2)))
		initial_condition[int(Nc) + 1] = 1.
		# print(initial_condition)

		ReL = np.array(simulation['linblad_real_part'])
		ImL = np.array(simulation['linblad_imaginary_part'])
		L = ReL + 1j * ImL


		# for k in range(int(Nc)):
		# 	for l in range(int(Nc)):
		# 		print(sum([L[200][j,j,k,l] for j in range(int(Nc))]))


		# print(ReL[0])
		pool = multiprocessing.Pool(processes = meta.CPU_COUNT)
		L_superoperator = pool.map(meta_functions.superoperator,
												tqdm.tqdm(
														zip(
															L, 
															[int(Nc)]*len(L)
															)
														)
													)
		# print(L[0])
		# Bloch-Redfield master equation is solved here:
		solution = master_eq_solve(initial_condition, L_superoperator, tvals, 
																		int(Nc))
		# print(tvals)
		# unpack components of the solution into the JSON dict
		simulation['rho_real_part'] = np.array(
											[np.real(p) for p in solution]).tolist()	
		simulation['rho_imaginary_part'] = np.array(
											[np.imag(p) for p in solution]).tolist()	



		# put the dict back into the JSON file
		print("uploading data into the JSON file...")

		start = time()
		with open(simulation_filename, 'w') as file:
			json.dump(simulation, file)
		end = time()
		print("Process complete. uploaded the solution of the master equation into the JSON file in {} seconds.".format(end-start))




