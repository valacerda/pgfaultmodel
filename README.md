# MMC-HVDC PG fault analytical expressions

The codes used to obtain the expressions in the paper:
_Approximated Analytical Model of Pole-to-ground Faults in Symmetrical Monopole MMC-HVDC Systems_ are provided.

## Usage

All codes were written in Matlab syntax.
The codes can be used separately.

Use **Kirchhoff_to_ODEs.m** to start with the system ODEs

Use **ODE_numerically.m** to solve the system of ODEs numerically

Use **ODEs_to_Laplace.m** to follow the derivation of the system solution in the frequency domain.

Use **Laplace_to_timedomain.m** to follow the simplifications in frequency domain and the resulting expressions in time domain

Use **resonant_frequencies.m** to follow the derivation of the system's undamped resonant frequencies
