# from sys import argv
# from meta_functions import map_level

import multiprocessing
import tqdm
import numpy as np
from scipy import linalg
import integrate
import matplotlib.pyplot as plt 
import meta
import meta_functions
from functools import partial


# tQA, I, J, K, N, Nc, step, window_size, num_samples, decoherence, CPU = [5*10**(-9), 0.2, 0.3, 1, 6, 6, 0.001, 10, 10000, 0, "Home"]
# args = [tQA, I, J, K, N, Nc, step, window_size, num_samples, decoherence, CPU]
# for j in range(len(args)-1):
#     args[j] = float(args[j])

# decoherence = int(decoherence)
# # args for constructing F-model Hamiltonian
# H_args = [I, J, K, int(N)]

# # args pertaining to numerical implementation of the QME
# N_args = [int(Nc), [step, window_size, int(num_samples)]]
# # print(int(num_samples))
# svals, bad_svals, tvals, LZ_probability = meta_functions.generate_discretization(
# 		tQA, H_args, N_args[1])

# hamiltonians = [meta_functions.generate_hamiltonian(svals[j], H_args) 
# 					for j in range(5960, 5966)]

# Nc = 6

# eigendata = [linalg.eigh(h, eigvals = (0, Nc - 1)) for h in hamiltonians]

# # print eigendata[0]
# ekets = [data[1] for data in eigendata]

# for j in range(len(ekets)):
# 	print ("ekets[{}][:,1] is:".format(j))
# 	print ekets[j][:,1]
# 	print ("ekets[{}][:,2] is:".format(j))
# 	print ekets[j][:,2]

# def foo(args):
# 	x, y = args
# 	x = float(x)
# 	args = [x, y]
# 	return x + y

# print(foo([1,2]))

# dt = 0.01
# tvals = np.arange(0, 1, dt)
# initial_condition = [1,1]
# k = -100

# L_superoperator = [k * np.identity(2, dtype = complex)]*len(tvals)
# rho = integrate.master_eq_solve(initial_condition, L_superoperator, tvals)

# plt.plot(tvals, list(map(np.exp, k*tvals)), tvals[:-2], rho[:-2,0])
# plt.show()


# for j in range(len(tvals)):
# 	if j % 2 == 1:
# 		list_of_t[j] = list_of_t[j - 1] + 2 * gamma * dt

# print (len(list_of_t))
# print (len(list_of_t)/25)
# print(list_of_t)




# print("hello.")
# print(list_of_linblads[0])

# list_of_rho_Explicit_Euler = ODE_integrator.run_time_evolution(rho_0, list_of_linblads, list_of_t, "Explicit Euler")
# list_of_rho_TR = ODE_integrator.run_time_evolution(rho_0, list_of_linblads, list_of_t, "TR")
# list_of_rho_BDF2 = ODE_integrator.run_time_evolution(rho_0, list_of_linblads, list_of_t, "BDF2")


# pylab.plot(list_of_t, list(map(np.exp, k*list_of_t)), '-k',label = 'Exact Solution')
# pylab.plot(list_of_t[:-2], list_of_rho_TRBDF2[:-2,0], '-m',label = 'TR-BDF2')
# pylab.plot(list_of_t[:-2], list_of_rho_TR[:-2,0], '-b',label = 'TR')
# pylab.plot(list_of_t[:-2], list_of_rho_BDF2[:-2,0], '-r',label = 'BDF2')


# pylab.legend(loc='upper left')


# pylab.show()
# terminal_input = argv[1:]
# print(terminal_input)
# args = map_level(float, terminal_input, 2)

# args = map(float, sys.argv[1:])
# print(args)
# import time

# N = 2**6
# Nc = 5
# H = np.random.rand(N, N)

# start = time.time()
# np_evals, np_ekets = np.linalg.eigh(H)
# end = time.time()
# print("numpy solver took {} seconds.".format(end-start))

# start = time.time()
# sp_evals, sp_ekets = linalg.eigh(H, eigvals = (0, Nc))
# end = time.time()
# print("scipy solver took {} seconds.".format(end-start))

N = 2**6
Nc = 5
H = np.random.rand(N, N)

# def diagonalization_func(Nc):
# 	return partial(linalg.eigh, eigvals = (0, Nc - 1))

evals, ekets = meta_functions.diagonalization_func(Nc)(H)

print(evals[0])


