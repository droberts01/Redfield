#---------------
# meta_functions.py 
# Provides global functions for meta, generate, and redfield.py
import meta
import numpy as np
from bisect import bisect
from scipy import linalg
from functools import partial

# import matplotlib.pyplot as plt

# Kronecker delta
def delta(i, j):
	if i == j: 
		return 1
	else:
		return 0

# Implements tensor product of a list of matrices
def tensor(matrices):
    return reduce(lambda x, y: np.kron(x,y), matrices)

# Define generators of the Pauli matrix algebra for N qubits
def X(q, N):
    if q >= N:
        print ("Error. Pauli matrix over-indexed")
    else:
        def matrices(q, k):
            if k == q:
                return meta.sigma_x
            else: 
                return np.identity(2)
        return tensor([matrices(q, k) for k in range(N)])

def Y(q, N):
    if q >= N:
        print ("Error. Pauli matrix over-indexed")
    else:
        def matrices(q, k):
            if k == q:
                return meta.sigma_y
            else: 
                return np.identity(2)
        return tensor([matrices(q, k) for k in range(N)])


def Z(q, N):
    if q >= N:
        print ("Error. Pauli matrix over-indexed")
    else:
        def matrices(q, k):
            if k == q:
                return meta.sigma_z
            else: 
                return np.identity(2)
        return tensor([matrices(q, k) for k in range(N)])


# Defines F-model Hamiltonian
def generate_hamiltonian(args):
	sval, H_args = args
	I, J, K, N = H_args

	def coupling(q):
		if q >= N - 1:
			print ("ERR: coupler is over-indexed")
		else:
			if N % 2 == 0:
				if q in [N/2 - 1, N/2]:
					return -J
				else:
					return -K
			else:
				if q in [int(N/2) - 1, int(N/2)]:
					return -J
				else:
					return -K
        

    # Defines the computational primitive; a 1D spin glass
	qs = range(N)
	bulk_qs = range(N-1)

	def bulk_term(q): 
		return coupling(q) * Z(q, N)* Z(q+1, N)

	boundary_term = I * Z(N-1, N) * Z(0, N)

	H_P = (sum([bulk_term(q) for q in bulk_qs]) + boundary_term)
	H_D = sum([X(q, N) for q in qs])

	return 2 * np.pi * (meta.A(sval)/2 * H_D + meta.B(sval)/2 * H_P)



def get_index(ordered_list, val):
    test = 0
    while val < ordered_list[test]:
        test += 1
    return test

# Generates time steps for the quantum master equation in
# a manner that dynamically increases resolution at any level crossing.
def generate_discretization(tQA, H_args, N_args_1):
	I, J, K, N = H_args
	step, window_size, num_samples = N_args_1

	svals = np.arange(0, 1, step)
	tvals = [tQA * s for s in svals]
	bad_svals = [0]
	# Test if bottleneck exists
	if I * K > J**2:

		# Solve for bottleneck location
		num = (K**2 - J**2)*(J**2 - I**2)
		den = (I**2 + K**2 - 2*J**2)
		gamma_b = (1./I) * num / den
		gamma_vals = map(lambda s: meta.A(s)/meta.B(s), svals)
		b_index = get_index(gamma_vals, gamma_b)
		s_b = svals[b_index]

		# Solve for minimum gap
		num = J**2 - I**2
		den = K**2 - J**2
		r_b = (K/I) * num / den
		min_gap = 0.05 * (r_b)**((N - 5.)/ 2.)


		# Generate LZ probability
		H_b = generate_hamiltonian([s_b, H_args])
		evals_b, ekets_b = linalg.eigh(H_b, eigvals = (0, 3))
		wmin = evals_b[3] - evals_b[0]
		H_b1, H_b2 = [generate_hamiltonian([s_b - 5 * min_gap, H_args]),
						generate_hamiltonian([s_b - 4 * min_gap, H_args])]
		evals_b1, ekets_b1 = linalg.eigh(H_b1, eigvals = (0, 3))
		evals_b2, ekets_b2 = linalg.eigh(H_b2, eigvals = (0, 3))
		wdot = np.abs((evals_b1[3] - evals_b1[0]) - (evals_b2[3] - evals_b2[0]))\
					/ min_gap
		# LZ_arg = -np.pi * wmin**2 * tQA / (2 * wdot)
		# LZ_probability = 1 - np.exp(LZ_arg)


		# Sample points according to bottleneck size and location
		middle_min = max(s_b - window_size * min_gap, 0)
		middle_max = min(s_b + window_size * min_gap, 1)
		middle_svals = list(np.arange(middle_min, middle_max, 
								(middle_max - middle_min)/float(num_samples)))
		
		start_min = 0
		start_max = middle_min
		start_svals = list(np.arange(start_min, start_max, step))


		end_min = middle_max
		end_max = 1
		end_svals = list(np.arange(end_min, end_max, step))


		# if middle_svals does not reach the boundaries of the
		# annealing, sample the boundaries at ordinary resolution.		
		start_used = (middle_min > 0)
		end_used = (middle_max < 1)
		if start_used:
			svals = start_svals + middle_svals
			if end_used:
				svals += end_svals
		elif not start_used:
			if end_used:
				svals = middle_svals + end_svals
		
	# Remove duplicates
	svals = sorted(list(set(svals)))
	svals = np.array(svals)

	# Store regions where the resolution in svals jumps sharply
	# in bad_svals
	for index in range(1,len(svals)-1):
		forward_step = round(svals[index + 1] - svals[index], 
												meta.STEP_PRECISION)
		backward_step = round(svals[index] - svals[index - 1], 
												meta.STEP_PRECISION)
		if forward_step != backward_step:
			bad_svals += [index]

	bad_svals += [len(svals) - 1]

	# Remove duplicates
	bad_svals = sorted(list(set(bad_svals)))
	tvals = [tQA * s for s in svals]
	# print([svals[b-1:b+2] for b in bad_svals[1:-1]])

	# print(svals)
	return svals, bad_svals, tvals


# Noise spectral density of the qubits in the 
# D Wave 2X quantum annealer
def td_bath_coupling(sval):
	scaling_factor = meta.B(sval)/meta.B(1)
	return scaling_factor * meta.BATH_COUPLING


def S(sval):
	def S_func(w):
		if abs(w) > meta.BATH_CUTOFF_FREQ:
			return 0.0
		elif w != 0:
			num = td_bath_coupling(sval) * w
			den = 1 - np.exp(-meta.BETA * w)
			return num / den
		else:
			return td_bath_coupling(sval) / meta.BETA
	return S_func


def superoperator(args):
	rk4_tensor, Nc = args
	def components(n, m):
		i, j = divmod(n, Nc)
		k, l = divmod(m, Nc)
		return rk4_tensor[i, j, k, l]

	# test = np.array([[components(n, m) for m in range(Nc**2)]
	# 									for n in range(Nc**2)])
	# print(test.shape)
	return np.array([[components(n, m) for m in range(Nc**2)]
										for n in range(Nc**2)])

def prod(iterable):
    return reduce((lambda x, y: x*y), iterable)

def diagonalization_func(Nc):
	return partial(linalg.eigh, eigvals = (0, Nc - 1))

def matrix(args):
	vector, Nc = args
	def components(n, m):
		i = n * Nc + m
		return vector[i]

	return np.array([[components(n, m) for m in range(Nc)]
										for n in range(Nc)])

# Tests of generate_hamiltonian()
# print(generate_hamiltonian(.32, [.2, .3, 1, 5])[:4,:4])

# Tests of generate_discretization()
# print(generate_discretization(5*10**(-6), [.2, .23, 1, 5], [.1, .2, 20])[0][24:27])
# print(generate_discretization(5*10**(-6), [.2, .3, 1, 8], [.01, 5, 20000])[1])

# Test of S()
# S_func = S(1)
# print(S_func(.1))
# x = np.arange(-10**6, 10**6, 10**3)
# y = map(S_func, x)
# plt.plot(x, y)
# plt.show()

