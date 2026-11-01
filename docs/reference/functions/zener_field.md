# `imr.physics.zener_field`

`imr.physics.zener_field` returns the time derivative of the local polymer stress difference `sp = tau_p,rr - tau_p,thetatheta` of the linear Maxwell, Jeffreys, Zener, 2nd-order Fung Zener and Mooney-Rivlin Zener models (stress 5–8). Both solvers carry `sp` as a field when medium heat transfer is on, because the heating (`imr.physics.stress_power`) needs the local stress, not only its integral `Z1`.

For the algorithm and assumptions, see the [topic guide](../../theory/stress-heating.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
