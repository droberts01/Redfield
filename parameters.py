"""
Original Developer: David Roberts
Purpose of Module: To specify input parameters for the time-dependent Redfield solver
Last Modified: 5/30/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
import pandas as pd


# Specifies mode of operation
# MODE = LANL
MODE = "GOOGLE"


# Fundamental Constants in SI Units
KB = 1.38065 *(10**(-23))
HBAR = 6.62607*(10**(-34))

# Operating Parameters
if MODE == "GOOGLE":
	ANNEALING_SCHEDULE = pd.read_csv('~/Documents/LANLA/DWaveAnnealingSchedule_Google.csv', sep=',',header=None)
elif MODE == "LANL":
	ANNEALING_SCHEDULE = pd.read_csv('~/Documents/LANLA/DWaveAnnealingSchedule.csv', sep=',',header=None)
else:
	print "ERR"

ANNEALING_PARAMETER = ANNEALING_SCHEDULE[0]
DRIVER_COEFFICIENT = [(10**9)*coefficient for coefficient in ANNEALING_SCHEDULE[1]]
PROBLEM_COEFFICIENT = [(10**9)*coefficient for coefficient in ANNEALING_SCHEDULE[2]]
ETA_MRT = 0.24
BATH_CUTOFF_TIME = 10**(-40)
TIME_STEP = ANNEALING_TIME*(ANNEALING_PARAMETER[1]-ANNEALING_PARAMETER[0])
INITIAL_DENSITY_MATRIX = [[1,0],[0,0]]


def ETA(s):
	def B(s):
		return PROBLEM_COEFFICIENT[s]
    return ETA_MRT*(B[s]/B[-1])

def SPECTRAL_DENSITY(s, frequency): 
    numerator =  (HBAR**2)*ETA(s)*frequency*np.exp(-np.abs(frequency)*BATH_CUTOFF_TIME)
    denominator = 1-np.exp(-(HBAR*frequency)/(KB*OPERATING_TEMPERATURE))
    return numerator/denominator



if MODE == "GOOGLE":
# Parameters of Google Computational Primitive
	h1 = .4
	h2 = -1
	J = 1
	NUM_QUBITS = 16
	ANNEALING_TIME = 5*(10**(-6)))
	OPERATING_TEMPERATURE = 15.5*(10**(-3))
	NUM_STATES = 2**NUM_QUBITS
	NUM_STATES_CUTOFF = 2
	cluster1 = []
	for i in range(4):
		for j in range(4)
			cluster1 = cluster1 + [[0 + i, 4 + j]]
	cluster2 = [[pair[0]+8,pair[1]+8] for pair in cluster1]
	pairs_connecting_cluster1_and_cluster_2 = [[4 + j, 12 + j] for j in range(4)]
	PAIR_SET = cluster1 + cluster2 + pairs_connecting_cluster1_and_cluster_2


elif MODE == "LANL":
# Parameters of Computational Primitive (F-model) 
	I = .2
	J = .3
	K = 1.0
	NUM_QUBITS = 5
	ANNEALING_TIME = 5*(10**(-6)))
	OPERATING_TEMPERATURE = 15.5*(10**(-3))
	NUM_STATES = 2**NUM_QUBITS
	NUM_STATES_CUTOFF = 2
else:
	print "ERR"







