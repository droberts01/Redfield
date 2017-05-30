# REDFIELD.py: Time-Dependent Redfield Solver for D Wave Simulations

from qutip import *
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy.linalg as linalg

# Parameters specific to the D Wave 2X onsite.
ANNEALING_SCHEDULE = pd.read_csv('~/Documents/LANLA/DWaveAnnealingSchedule.csv', sep=',',header=None)
ANNEALING_PARAMETER = ANNEALING_SCHEDULE[0]
DRIVER_COEFFICIENT = [(10**9)*x for x in ANNEALING_SCHEDULE[1]]
PROBLEM_COEFFICIENT = [(10**9)*x for x in ANNEALING_SCHEDULE[2]]

ETA_MRT = 0.24

# Fundamental Constants in SI Units
KB = 1.38065 *(10**(-23))
T = 15.5*(10**(-3))
HBAR = 6.62607*(10**(-34))

#Define Pauli matrices. e.g. X(actingQubit) is the bit flip operator on actingQubit.

def X(acting_qubit, num_qubits):
    qubits = range(num_qubits)
    if acting_qubit >= num_qubits:
        print "Error. Pauli matrix over-indexed"
    else:
        def X_tensor(acting_qubit, qubit):
            if qubit == acting_qubit:
                return sigmax()
            else: 
                return identity(2)
        return tensor([X_tensor(acting_qubit, qubit) for qubit in qubits])

def Y(acting_qubit, num_qubits):
    qubits = range(num_qubits)
    if acting_qubit >= num_qubits:
        print "Error. Pauli matrix over-indexed"
    else:
        def Y_tensor(acting_qubit, qubit):
            if qubit == acting_qubit:
                return sigmay()
            else: 
                return identity(2)
        return tensor([Y_tensor(acting_qubit, qubit) for qubit in qubits])


def Z(acting_qubit, num_qubits):
    qubits = range(num_qubits)
    if acting_qubit >= num_qubits:
        print "Error. Pauli matrix over-indexed"
    else:
        def Z_tensor(acting_qubit, qubit):
            if qubit == acting_qubit:
                return sigmaz()
            else: 
                return identity(2)
        return tensor([Z_tensor(acting_qubit, qubit) for qubit in qubits])


#Define D Wave Hamiltonian
def bulk_coupling(J, K, qubit_index, num_qubits):
    if qubit_index >= num_qubits - 1:
        print "Error. coupler is over-indexed"
    else:
        if num_qubits % 2 == 0:
            if qubit_index in [num_qubits/2 - 1, num_qubits/2]:
                return -J
            else:
                return -K
        else:
            if qubit_index in [num_qubits/2 - 1, num_qubits/2]:
                return -J
            else:
                return -K
               
def d_wave_hamiltonian(s, I, J, K, num_qubits):
    qubits = range(num_qubits)
    s_rescaled = int(max(0, s*round(len(ANNEALING_PARAMETER)) - 1))
    bulk_terms = sum([bulk_coupling(J, K, qubit, num_qubits)*Z(qubit, num_qubits)*Z(qubit + 1, num_qubits) for qubit in range(num_qubits - 1)])
    boundary_term = I*Z(qubits[-1], num_qubits)*Z(qubits[0], num_qubits)
    raw_problem_hamiltonian = bulk_terms + boundary_term
    raw_driver_hamiltonian = sum([X(qubit, num_qubits) for qubit in qubits])
    problem_hamiltonian = 2*np.pi*(PROBLEM_COEFFICIENT[s_rescaled]/2)*raw_problem_hamiltonian
    driver_hamiltonian = 2*np.pi*(DRIVER_COEFFICIENT[s_rescaled]/2)*raw_driver_hamiltonian
    return problem_hamiltonian + driver_hamiltonian
                    
def eta(s):
    s_rescaled = int(max(0, s*round(len(ANNEALING_PARAMETER)) - 1))
    return ETA_MRT*(B[s_rescaled]/B[-1])


# Spectral Density

def S(s, omega): 
    numerator =  (hbar**2)*eta(s)*omega*np.exp(-np.abs(omega)*10**(-40))
    denominator = 1-np.exp(-(HBAR*omega)/(KB*T))
    return numerator/denominator

def delta(i,j):
    if i == j:
        return 1
    else:
        return 0

def O_tensor(O_matrix):
    states = range(len(O_matrix[0]))
    def O(s,i,j,k,l):
        return delta(j,l)*delta(i,k)*O_matrix[s][i][j]
    return [[[[[O(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]

def R_tensor(eigenvecs_Qobj, O_matrix, num_qubits):
    states = range(len(eigenvecs_Qobj[0]))
    qubits = range(num_qubits)
    def C(s,i,j,k,l):
        def Z_elt(s,mu,i,j):
            Z(mu,num_qubits).matrix_element(eigenvecs_Qobj[s][i],eigenvecs_Qobj[s][j])
        return sum([Z_elt(s,mu,i,k)*Z_elt(s,mu,j,l) for mu in qubits])
    
    C_tensor = [[[[[C(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]
    
    def gamma_plus(s,i,j,k,l):
        return S(-O_matrix[s][k][l])*C_tensor[s][i][j][k][l]
    def gamma_minus(s,i,j,k,l):
        return S(O_matrix[s][i][j])*C_tensor[s][i][j][k][l]

    def R1(s,i,j,k,l):
        return sum([delta(j,l)*gamma_plus(s,i,k,n,n) for n in range(2**num_qubits)])

    def R2(s,i,j,k,l):
        return sum([delta(j,l)*gamma_plus(s,i,k,n,n) for n in range(2**num_qubits)])

    def R3(s,i,j,k,l):
        return gamma_plus(s,i,j,k,l)

    def R4(s,i,j,k,l):
        return gamma_minus(s,i,j,k,l)
    
    R1_tensor = [[[[[R1(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]
    R2_tensor = [[[[[R2(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]
    R3_tensor = [[[[[R3(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]
    R4_tensor = [[[[[R4(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]
    
    return R1_tensor + R2_tensor - R3_tensor - R4_tensor


def M_tensor(eigenvecs_Qobj, num_qubits, ANNEALING_TIME):
    states = range(len(eigenvecs_Qobj[0]))
    TIME_STEP = ANNEALING_TIME*(ANNEALING_PARAMETER[1]-ANNEALING_PARAMETER[0])
    def basis_dragging_term(s,i,j):
        bra = eigenvecs_Qobj[s + 1][i] + eigenvecs_Qobj[s - 1][i]
        ket = eigenvecs_Qobj[s + 1][i] - eigenvecs_Qobj[s - 1][i]
        return identity(2**num_qubits).matrix_element(bra, ket)/(4*TIME_STEP)
    
    def M(s,i,j,k,l):
        if s == 0 or s == len(ANNEALING_PARAMETER)-1:
            return 0
        else:
            return delta(i,k)*basis_dragging_term(s,l,j) + delta(j,l)*basis_dragging_term(s,i,k)
    return [[[[[M(s,i,j,k,l) for l in states] for k in states] for j in states] for i in states] for s in ANNEALING_PARAMETER]


def linbladian(I,J,K, eigenvecs_Qobj, num_qubits, ANNEALING_TIME, STATE_CUT_OFF):
    hamiltonian = [d_wave_hamiltonian(s, I, J, K, num_qubits) for s in ANNEALING_PARAMETER]
    eigenvalues = hamiltonian[s].eigenenergies()
    eigenvecs = hamiltonian[s].eigenstates()
    idx = eigenvalues.argsort()[::-1]   
    eigenvalues_sorted = eigenvalues[idx]
    eigenvecs_sorted = eigenvecs[:,idx]
    eigenvalues_sorted_truncated = eigenvalues_sorted[::STATE_CUT_OFF]
    eigenvecs_sorted_truncated = eigenvecs_sorted[::STATE_CUT_OFF]
    
    O_matrix = [[eigenvalues_sorted[i]-eigenvalues_sorted[j] for i in range(len(eigenvalues))] for j in range(len(eigenvalues))]
    O_matrix_truncated = O_matrix[::STATE_CUT_OFF,::STATE_CUT_OFF]
    
    L1_tensor = -1j*O_tensor(O_matrix_truncated)
    L2_tensor = R_tensor(eigenvecs_sorted, O_matrix, num_qubits)
    L3_tensor = M_tensor(eigenvecs_sorted_truncated,num_qubits, ANNEALING_TIME)
    
    return L1_tensor - (L2_tensor + L3_tensor)
    
# 