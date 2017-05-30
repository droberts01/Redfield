"""
Original Developer: David Roberts
Purpose of Module: outputs Redfield tensor
Last Modified: 5/30/17
Last Modified By: David Roberts
Last Modification Purpose: fixed function naming
"""


# Standard Modules:
from qutip import *
import csv
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as linalg

# Custom Modules:
import parameters


# Spectral Density

def S(s, omega): 
    numerator =  (hbar**2)*eta(s)*omega*np.exp(-np.abs(omega)*parameters.BATH_CUTOFF_TIME)
    denominator = 1-np.exp(-(parameters.HBAR*omega)/(parameters.KB*parameters.OPERATING_TEMPERATURE))
    return numerator/denominator


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