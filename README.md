# MMC-HVDC PG fault analytical expressions

This repository provides the codes used to obtain the expressions in the paper: _V. A. Lacerda, R. M. Monaro, D. Campos-Gaona, R. Pena-Alzola, D. V. Coury "An Approximated Analytical Model of Pole-to-ground Faults in Symmetrical Monopole MMC-HVDC Systems. IEEE Journal on Emerging and Selected Topics in Power Electronics, 2020._"

## Usage

All codes were written in Matlab syntax.
The codes can be used separately.

Use **Kirchhoff_to_ODEs.m** to start with the system ODEs

Use **ODE_numerically.m** to solve the system of ODEs numerically

Use **ODEs_to_Laplace.m** to follow the derivation of the system solution in the frequency domain.

Use **Laplace_to_timedomain.m** to follow the simplifications in frequency domain and the resulting expressions in time domain

Use **resonant_frequencies.m** to follow the derivation of the system's undamped resonant frequencies
