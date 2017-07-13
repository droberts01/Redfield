"""
Original Developer: David Roberts
Purpose of Module: provides numerical methods for generating the integral
solution to the ODE as performed in td_linblad_solver.py.
Last Modified: 6/20/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""


import numpy as np
import parameters


GAMMA = 2. - np.sqrt(2)
NUM_STATES_CUTOFF = parameters.NUM_STATES_CUTOFF

# 1st-Order Adams-Bashforth (Explicit Euler) Method

def Explicit_Euler_update(rho_n, L_n, dt_n):
	return rho_n + dt_n * np.matmul(L_n, rho_n)

# 2nd-Order TR-BDF2 (Backwards Differentiation) Method

# def TR_BDF2_update(rho_n, L_n, L_n1, dt_n):
# 	I = np.identity(len(rho_n), dtype = complex)
# 	B = np.linalg.inv(I - GAMMA * dt_n * L_n / 2.)
# 	B2 = np.linalg.inv(I + GAMMA * dt_n * L_n / 2.)
# 	return np.matmul(B, 
# 		np.matmul(
# 			(1./(GAMMA*(2-GAMMA))) * np.matmul(B, B2) - ((1 - GAMMA)**2)/(GAMMA * (2 - GAMMA)) * I, 
# 					rho_n)  
# 		)


def Implicit_Euler_update(rho_n, L_n1, dt_n):
	I = np.identity(len(rho_n), dtype = complex)
	B = np.linalg.inv(I - dt_n * L_n1)
	return np.matmul(B, rho_n)


def TR_update(rho_n, L_n, L_n1, dt_n):
	I = np.identity(len(rho_n), dtype = complex)
	B1 = np.linalg.inv((I - (dt_n/2.)*L_n1))
	B2 = (I + (dt_n/2.)*L_n)
	return np.matmul(B1, np.matmul(B2, rho_n))

def BDF2_update(rho_n, rho_n_minus_1, L_n1, dt_n):
	I = np.identity(len(rho_n), dtype = complex)
	B = np.linalg.inv(1.5*I - dt_n*L_n1)
	return np.matmul(B, 2* rho_n - .5* rho_n_minus_1)

def run_time_evolution(rho_0, list_of_linblads, list_of_t, method):
	rho = np.array([rho_0]*len(list_of_t), dtype=complex)
	time_indices = range(len(list_of_t))
	for time_index in time_indices[:-2]:
		rho_n = rho[time_index]
		rho_n_minus_1 = rho[time_index - 1]
		L_n = list_of_linblads[time_index]
		L_n1 = list_of_linblads[time_index+1]
		dt_n = list_of_t[time_index + 1] - list_of_t[time_index]
		if method == "Explicit Euler":
			rho[time_index + 1] = Explicit_Euler_update(rho_n, L_n, dt_n) 
		elif method == "TR-BDF2":
			if time_index%2 ==0:
				rho[time_index + 1] = TR_update(rho_n, L_n, L_n1, dt_n)
			else:
				rho[time_index + 1] = BDF2_update(rho_n, rho_n_minus_1, L_n1, dt_n)
		elif method == "Implicit Euler":
			rho[time_index + 1] = Implicit_Euler_update(rho_n, L_n1, dt_n) 
		elif method == "TR":
			rho[time_index + 1] = TR_update(rho_n, L_n, L_n1, dt_n)
		elif method == "BDF2":
			rho[time_index + 1] = BDF2_update(rho_n, rho_n_minus_1, L_n1, dt_n)
		else:
			print("ERR. Method not specified.")

		if time_index%(len(list_of_t)/20) == 0:
			print("time_index is {}".format(time_index))
			print ("time_step is {}".format(dt_n))
			print("tr(rho[time_index]) is")
			print(np.transpose(rho[time_index-1,np.newaxis]))
			print(sum([rho[time_index-1,(NUM_STATES_CUTOFF+1)*j]  for j in range(NUM_STATES_CUTOFF)]))

	return rho
				
