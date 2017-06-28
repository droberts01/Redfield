"""
Original Developer: David Roberts
Purpose of Module: exports linblad.csv according to Redfield approximation
Last Modified: 6/5/17
Last Modified By: David Roberts
Last Modification Purpose: Created module
"""

# Standard Modules:
from qutip import *
import numpy as np
from multiprocessing import Pool
import tqdm
from pathlib import Path

# Custom Modules:
import parameters
import redfield
import primitive
import helper
import time

# Load parameters
list_of_s = parameters.S_VALUES
list_of_t = parameters.LIST_OF_TIMES





globalstart = time.time()

root_filepath = parameters.ROOT_FILEPATH
filename = parameters.FILENAME

print("testing to see if linblad data already exists...")
my_file_reals = Path(root_filepath+filename+'reals_v2.csv')
my_file_imags = Path(root_filepath+filename+'imags_v2.csv')

# WARNING: DEBUGGING PURPOSES ONLY
# if True:
if my_file_reals.is_file() and my_file_imags.is_file():
	print ("linblad for {} settings exists. To re-run, comment-out the conditional in linblad_v2.py.".format(filename))

else:
	print ("linblad for {} settings does not exist yet. Producing new linblads...".format(filename))
		


	# Program runs here:
	print ("computing hamiltonians...")
	start = time.time() 
	pool = Pool(processes=8)
	list_of_hamiltonians_Qobj = pool.map(primitive.hamiltonian, tqdm.tqdm(list_of_s))
	list_of_hamiltonians = [hamiltonian.full() for hamiltonian in list_of_hamiltonians_Qobj]
	end = time.time()
	print ("computed hamiltonians in {} seconds.".format(end-start))


	print ("computing redfield tensors...")
	start = time.time() 
	# list_of_redfield_output = map(redfield.compute_redfield_tensor, zip(list_of_s, list_of_hamiltonians_Qobj))
	pool = Pool(processes=8)
	list_of_redfield_output = pool.map(redfield.compute_redfield_tensor_v2, tqdm.tqdm(zip(list_of_s, list_of_hamiltonians)))

	list_of_redfield_tensors, list_of_eigenstates = zip(*list_of_redfield_output)
	end = time.time()

	print ("computed redfield tensors in {} seconds.".format(end-start))


	print ("truncating eigenstates...")
	list_of_eigenstates = map(helper.truncate, list_of_eigenstates)
	start = time.time() 

	print ("computing diabatic tensors...")
	pool = Pool(processes=8)
	dt = np.diff(list_of_t)
	eigenstate_window = [ [[0],[0]] ]+ [[list_of_eigenstates[i-1],list_of_eigenstates[i+1]] 
										for i in range(1,len(list_of_eigenstates)-1)] +[ [[0],[0]] ]

	list_of_diabatic_tensors = pool.map(redfield.compute_diabatic_tensors, 
				tqdm.tqdm(zip(eigenstate_window, dt, range(len(list_of_t)))))



	pool = Pool(processes=8)
	list_of_linblads = pool.map(helper.get_sum_tensors, 
							zip(list_of_redfield_tensors, list_of_diabatic_tensors))

	end = time.time()
	print ("computed diabatic tensors in {} seconds.".format(end-start))


	dim = parameters.NUM_STATES_CUTOFF
	pool = Pool(processes=8)
	list_of_compact_linblads = pool.map(helper.get_compact_compact_tensor_from_tensor, 
										zip(list_of_linblads, [dim]*len(list_of_s)))

	# print(list_of_compact_linblads[0])
	# print(np.imag(list_of_compact_linblads[0]))

	# Export array as csv
	list_of_compact_linblads_reals = np.array(map(np.real, list_of_compact_linblads))
	list_of_compact_linblads_imags = np.array(map(np.imag, list_of_compact_linblads))

	# print(list_of_compact_linblads_imags[0])

	np.savetxt(root_filepath+filename+'reals_v2.csv', list_of_compact_linblads_reals, delimiter=",")
	np.savetxt(root_filepath+filename+'imags_v2.csv', list_of_compact_linblads_imags, delimiter=",")


globalend = time.time()
print ("linblad_v2.py complete. process took {} seconds.".format(globalend-globalstart))




