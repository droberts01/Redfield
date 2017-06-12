"""
Original Developer: David Roberts
Purpose of Module: To specify input parameters for the time-dependent Redfield solver
Last Modified: 6/5/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
import pandas as pd
import numpy as np

# Specifies mode of operation
# MODE = "LANL"
MODE = "NASA"


# Fundamental Constants in SI Units
KB = 1.38065 *(10**(-23))
HBAR = 6.62607*(10**(-34))

# Operating Parameters
if MODE == "GOOGLE":
	ANNEALING_SCHEDULE = pd.read_csv('~/Documents/LANLA/DWaveAnnealingSchedule_Google.csv', sep=',',header=None)
elif MODE == "NASA":
	ANNEALING_SCHEDULE = pd.read_csv('~/Documents/LANLA/DWaveAnnealingSchedule.csv', sep=',',header=None)
else:
	print "ERR"

ANNEALING_TIME = 5*(10**(-6))
ANNEALING_PARAMETER = ANNEALING_SCHEDULE[0]
DRIVER_COEFFICIENT = [(10**9)*coefficient for coefficient in ANNEALING_SCHEDULE[1]]
PROBLEM_COEFFICIENT = [(10**9)*coefficient for coefficient in ANNEALING_SCHEDULE[2]]
BATH_COUPLING_MRT = 0.24
BATH_CUTOFF_TIME = 10**(-12)
NUM_STATES_CUTOFF = 4
INITIAL_DENSITY_MATRIX = [[0]*NUM_STATES_CUTOFF]*NUM_STATES_CUTOFF
INITIAL_DENSITY_MATRIX[0][0] = 1


# Defines discretization of Linblad ODE
S_VALUES = np.arange(0,1,.5)
LIST_OF_TIMES = [ANNEALING_TIME * s for s in S_VALUES]


def A(s):
	s_int = max(int(round(len(ANNEALING_PARAMETER) * s)) - 1, 0)
	return DRIVER_COEFFICIENT[s_int]

def B(s):
	s_int = max(int(round(len(ANNEALING_PARAMETER) * s)) - 1, 0)
	return PROBLEM_COEFFICIENT[s_int]




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
	NUM_QUBITS = 10
	OPERATING_TEMPERATURE = 15.5*(10**(-3))
	NUM_STATES = 2**NUM_QUBITS
else:
	print "ERR"







