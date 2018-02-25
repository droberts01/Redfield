# Redfield
A portable time-dependent Redfield solver for high-performance simulations of QAA

Calculation of transition rates suggests that incoherent tunneling accelerates at bottlenecks of QAA, serving as a considerable computational resource. This Redfield solver facilitates a search for signatures of this hardware acceleration in the D Wave 2X adiabatic quantum computer at NASA Ames.

## Getting Started
Our program "Redfield" only runs on Python 2, not Python 3. So make sure that all commands of the type

    python FILENAME.py

Are running Python 2, not Python 3.
Finally, install scipy and numpy, if you haven't done so.

As a last step, please specify the folder that you want to save/retrieve data from by setting the appropriate value for the variable meta.SAVE_LOCATION (this variable is defined in line 9 of Redfield/src/meta.py):

    # Detect local directory
    ROOT = os.path.dirname(os.path.realpath(__file__))
    SAVE_LOCATION = 'something/something/something'

You'll want to change this to the save location of your choice. For example, on a Mac, the following code would suffice:

    SAVE_LOCATION = '/Users/USERNAME/Documents/PLACE_I_WANT_TO_STORE_MY_OUTPUT'

## Running D Wave Annealing Simulations

To run the Redfield solver, first run 

    python generate.py tQA I J K N Nc step window_size num_samples include_low_freq_noise include_decoherence store_linblads save_location
    
The options here are as follows:

- tQA denotes the total anneal time in seconds (e.g. 5E-6)
- I,J,K specify the coupling strengths 
- N specifies the total number of qubits.

- Nc specifies the dimensions of the density matrix, which is Nc x Nc. Equivalently, this is the truncated Hilbert space dimension.
- step specifies the timestep size (in units of the rescaled time variable 0 < s < 1) to use for the simulation.
- window_size and num_samples specify how finely the bottleneck region is sampled.

Finally, we have some Boolean options (True or False):
- include_low_freq_noise specifies whether to include a 1/f component in the noise spectral density.
- include_decoherence specifies whether to include the effects of the bath or to just simulate the coherent evolution (i.e. time-dependent closed-system dynamics)
- store_linblads specifies whether to save the time-dependent linblad superoperator to the output file.

And, the very last option:
- save_location specifies where to save the file (if you are a new user, select "Home").

Finally, once generate.py has run, run solve.py:

    python solve.py tQA I J K N Nc step window_size num_samples include_low_freq_noise include_decoherence store_linblads save_location

The options here have to be exactly the same as the ones used in the prior generate command, otherwise the terminal will throw an error.

# Example

Type
   
    python generate.py 5E-6 .2 0.24 1 7 6 0.001 10 1000 1 1 0 'Home'

to generate the Bloch-Redfield quantum master equation, then run


    python generate.py 5E-6 .2 0.24 1 7 6 0.001 10 1000 1 1 0 'Home'

To solve this master equation.
