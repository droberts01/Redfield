"""
Original Developer: David Roberts
Purpose of Module: outputs solution to the time-dependent Linblad equation.
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
import redfield
import primitive
import linblad_solver



def linblad_superoperator(density_matrix):
	return sum()