#---------------
# meta.py 
# Provides global variables for meta, generate, and redfield.py
import os
import multiprocessing
from scipy.interpolate import interp1d
import numpy as np

# Detect local directory
ROOT = os.path.dirname(os.path.realpath(__file__))
SAVE_LOCATION = '/Users/Droberts/Documents/Dissipative/Redfield_data/'

# Detect local parallelization resources
CPU_COUNT = multiprocessing.cpu_count()
STEP_PRECISION = 5

# Operating conditions in the D Wave 2X
T = 0.0155
KB = 1.38065 *(10**(-23))
HBAR = 1.0545718*(10**(-34))
BETA = HBAR / (KB * T)

ANNEALING_SCHEDULE = np.transpose(
	np.loadtxt('data/DWaveAnnealingSchedule.csv', delimiter=','))
ANNEALING_PARAMETER, DRIVER_AMPLITUDE, PROBLEM_AMPLITUDE  = ANNEALING_SCHEDULE

# File is implicitly in units of GHz
DRIVER_AMPLITUDE *= 10**9
PROBLEM_AMPLITUDE *= 10**9

# Linearly interpolate driver and problem amplitudes
A = interp1d(ANNEALING_PARAMETER, DRIVER_AMPLITUDE, kind = 'linear')
B = interp1d(ANNEALING_PARAMETER, PROBLEM_AMPLITUDE, kind = 'linear')

# Time-Dependent (TD) noise strength in the D Wave 2X
BATH_COUPLING = 0.24
BATH_CUTOFF_FREQ = 10**11

# # Initial qubit state in QAA
# INITIAL_STATE = 

# Pauli matrices
sigma_x = np.array([[0,1],[1,0]])
sigma_y = np.array([[0,1j],[-1j,0]])
sigma_z = np.array([[1,0],[0,-1]])
