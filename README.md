# Redfield
Time-dependent Redfield solver for simulating D Wave physics

We analytically investigate the noise-susceptibility of the quantum annealing (QA) algorithm used by the D-Wave quantum adiabatic computer. To do this, we solve a special class of models which contain frustrated bottlenecks of the kind expected to dominate the runtime of QA for problem sizes larger than  N ~ 10^3 qubits. Calculation of transition rates suggests that incoherent tunneling accelerates at these bottlenecks, thermalizing these bottlenecks in an intermediate-N ("soft") regime, dominating coherent (quantum) tunneling in the large-N limit, and forming a computational advantage. This Redfield solver spearheads an ongoing search for experimental signatures of this new phenomenon in the D Wave 2X at NASA Ames.
