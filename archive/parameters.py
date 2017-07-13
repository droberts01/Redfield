"""
Original Developer: David Roberts
Purpose of Module: To specify input parameters for the time-dependent Redfield solver
Last Modified: 6/5/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
import numpy as np
import matplotlib.pyplot as plt
import os


ROOT_FILEPATH = "/Users/Droberts/Dropbox/Redfield Annealing/data/"


# Specifies mode of operation
# MODE = "LANL"
MODE = "NASA"


# Fundamental Constants in SI Units
KB = 1.38065 *(10**(-23))
HBAR = 1.0545718*(10**(-34))

# Operating Parameters
if MODE == "GOOGLE":
	ANNEALING_SCHEDULE = np.transpose(np.loadtxt('data/DWaveAnnealingSchedule_Google.csv', delimiter=','))
elif MODE == "NASA":
	ANNEALING_SCHEDULE = np.transpose(np.loadtxt('data/DWaveAnnealingSchedule.csv', delimiter=','))
else:
	print "ERR"

ANNEALING_TIME = 5*(10**(-6))
ANNEALING_PARAMETER = ANNEALING_SCHEDULE[0]
DRIVER_COEFFICIENT = [(10**9)*coefficient for coefficient in ANNEALING_SCHEDULE[1]]
PROBLEM_COEFFICIENT = [(10**9)*coefficient for coefficient in ANNEALING_SCHEDULE[2]]
BATH_COUPLING_MRT = 0.24
BATH_CUTOFF_TIME = 10**(-12)
NUM_STATES_CUTOFF = 8
INITIAL_DENSITY_MATRIX = np.array([[0]*NUM_STATES_CUTOFF]*NUM_STATES_CUTOFF)
INITIAL_DENSITY_MATRIX[0][0] = 1


# Defines discretization of Linblad ODE
S_VALUES = np.arange(0,1, 2*(10**(-4)))
LIST_OF_TIMES = [ANNEALING_TIME * s for s in S_VALUES]
DT = LIST_OF_TIMES[1]-LIST_OF_TIMES[0]




def A(s):
	s_rescaled = len(ANNEALING_PARAMETER) * s
	s_low = max(int(np.floor(s_rescaled)) - 1, 0)
	s_high = max(int(np.ceil(s_rescaled)) - 1, 0)
	if s_low == s_high:
		return DRIVER_COEFFICIENT[s_low]
	elif s_high == len(ANNEALING_PARAMETER):
		return DRIVER_COEFFICIENT[-1]
	else:
		return DRIVER_COEFFICIENT[s_low] + (s_rescaled - s_low)/(s_high - s_low) * (DRIVER_COEFFICIENT[s_high] - DRIVER_COEFFICIENT[s_low])

def B(s):
	s_rescaled = len(ANNEALING_PARAMETER) * s
	s_low = max(int(np.floor(s_rescaled)) - 1, 0)
	s_high = max(int(np.ceil(s_rescaled)) - 1, 0)
	if s_low == s_high:
		return PROBLEM_COEFFICIENT[s_low]
	elif s_high == len(ANNEALING_PARAMETER):
		return PROBLEM_COEFFICIENT[-1]
	else:
		return PROBLEM_COEFFICIENT[s_low] + (s_rescaled - s_low)/(s_high - s_low) * (PROBLEM_COEFFICIENT[s_high] - PROBLEM_COEFFICIENT[s_low])


# x = ANNEALING_PARAMETER
# y1 = map(A, ANNEALING_PARAMETER)
# y2 = map(B, ANNEALING_PARAMETER)

# plt.plot(x,y1,x,y2)
# plt.grid(True)
# plt.show()


def BATH_COUPLING(s):
    return BATH_COUPLING_MRT*(B(s)/B(1))


if MODE == "GOOGLE":
# Parameters of Google Computational Primitive
	h1 = .4
	h2 = -1
	J = 1
	NUM_QUBITS = 16
	OPERATING_TEMPERATURE = 15.5*(10**(-3))
	NUM_STATES = 2**NUM_QUBITS
	cluster1 = []
	for i in range(4):
		for j in range(4):
			cluster1 = cluster1 + [[0 + i, 4 + j]]
	cluster2 = [[pair[0]+8,pair[1]+8] for pair in cluster1]
	pairs_connecting_cluster1_and_cluster_2 = [[4 + j, 12 + j] for j in range(4)]
	PAIR_SET = cluster1 + cluster2 + pairs_connecting_cluster1_and_cluster_2


elif MODE == "NASA":
# Parameters of Computational Primitive (F-model) 
	I = .2
	J = .3
	K = 1.0
	NUM_QUBITS = 6
	OPERATING_TEMPERATURE = 15.5*(10**(-3))
	NUM_STATES = 2**NUM_QUBITS

	T = OPERATING_TEMPERATURE
	tQA = ANNEALING_TIME
	N = NUM_QUBITS
	N_c = NUM_STATES_CUTOFF
	dt = TIME_STEP = (LIST_OF_TIMES[2] - LIST_OF_TIMES[0])/2
	# GAMMA = 2 - np.sqrt(2)

	# for j in range(len(LIST_OF_TIMES)):
	# 	if j % 2 == 1:
	# 		LIST_OF_TIMES[j] = LIST_OF_TIMES[j - 1] + 2 * GAMMA * DT


	FILENAME = 'T='+str(T)+'___tQA='+str(tQA)+'___N='+str(N)+'___[I,J,K]='+str([I,J,K])+'___N_c='+str(N_c)+'___dt='+str(dt)+'___'


else:
	print "ERR"







