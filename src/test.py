# from sys import argv
# from meta_functions import map_level

import multiprocessing
import tqdm
import numpy as np
from scipy import linalg

# terminal_input = argv[1:]
# print(terminal_input)
# args = map_level(float, terminal_input, 2)

# args = map(float, sys.argv[1:])
# print(args)
import time

N = 2**11
Nc = 5
H = np.random.rand(N, N)


start = time.time()
np_evals, np_ekets = np.linalg.eigh(H)
end = time.time()
print("numpy solver took {} seconds.".format(end-start))

start = time.time()
sp_evals, sp_ekets = linalg.eigh(H, eigvals = (0, Nc))
end = time.time()
print("scipy solver took {} seconds.".format(end-start))





