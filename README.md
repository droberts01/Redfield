# Redfield
A portable time-dependent Redfield solver for high-performance simulations of QAA

Calculation of transition rates suggests that incoherent tunneling accelerates at bottlenecks of QAA, serving as a considerable computational resource. This Redfield solver facilitates a search for signatures of this hardware acceleration in the D Wave 2X adiabatic quantum computer at NASA Ames.

## Using Redfield

To run the Redfield solver, first run 

    generate.py

to generate the Bloch-Redfield quantum master equation, then run

    solve.py

To solve this master equation.
