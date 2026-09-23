# `f_display.m` Documentation

## Overview
`f_display` prints a summary of the run to the command window before time integration: the radial model, the constitutive model, the viscosity model, the thermal, mass-transfer and vapor settings, the solution method, and the main dimensionless numbers. Both solvers call it.

## Inputs
`radial`, `bubtherm`, `medtherm`, `masstrans`, `stress`, `spectral`, `nu_model`, `eps3`, `Pv_star`, `Re8`, `De`, `Ca`, `LAM`, and `mode`, the solver name printed as the solution method.

## Notes
- The constitutive label depends on the numbers as well as the model: for stress 1 and 2 with `Ca = Inf` it prints "Newtonian fluid"; for stress 3 and 4 it prints linear Zener, Maxwell, Jeffreys or the Kelvin–Voigt series from `Ca` and `LAM`; for stress 5 it prints upper-convected Maxwell when `LAM = 0` and Oldroyd-B otherwise.
