import numpy as np
from qutip import *

import primitive
import parameters
import redfield
import primitive
import helper

states = range(parameters.NUM_STATES_CUTOFF)

test_s = .2

test_hamiltonian_Qobj = primitive.hamiltonian(test_s)
test_hamiltonian = test_hamiltonian_Qobj.full()

R, eigenstates = redfield.compute_redfield_tensor_v2([test_s, test_hamiltonian])

def R_trace(k, l):
	return sum([R[i,i,k,l] for i in states])

print(np.array([[R_trace(i,j) for i in states] for j in states]))