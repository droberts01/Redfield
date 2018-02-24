# Redfield
A portable time-dependent Redfield solver for high-performance simulations of QAA

Calculation of transition rates suggests that incoherent tunneling accelerates at bottlenecks of QAA, serving as a considerable computational resource. This Redfield solver facilitates a search for signatures of this hardware acceleration in the D Wave 2X adiabatic quantum computer at NASA Ames.

## Getting Started
Our program "Redfield" only runs on Python 2, not Python 3. So make sure that all commands of the type

    python FILENAME.py

Are running Python 2, not Python 3.
Finally, install scipy and numpy, if you haven't done so.

## Running D Wave Annealing Simulations

To run the Redfield solver, first run 

    python generate.py

to generate the Bloch-Redfield quantum master equation, then run

    python solve.py

To solve this master equation.
