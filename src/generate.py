#---------------
# generate.py 
# Generates JSON dict with linblad for the Redfield quantum master equation
from sys import argv
import redfield

# Generate empty JSON dict with linblad
terminal_input = argv[1:]
args = [float(num) for num in terminal_input]


# Generate linblad for the Redfield quantum master equation
from time import time

start = time()
redfield.generate_json(args)
end = time()

print("Overall JSON runtime was {} seconds.".format(end-start))


