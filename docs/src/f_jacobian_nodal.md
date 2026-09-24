# `f_jacobian_nodal.m` Documentation

## Overview
`f_jacobian_nodal(t, X, jp)` returns the analytic Jacobian of the nodal right-hand side, for either differentiation scheme, for `bubtherm = 1`, with or without medium heat transfer (`medtherm`) and vapor mass transfer (`masstrans`). `f_imr_nodal` passes it to the stiff integrators when `jacobian = 2`.

## Method
- The field equations are differentiated exactly with the value/derivative pairs of `f_dual_ops` (shared with `f_jacobian_spectral`). Every intermediate quantity (temperature, vapor concentration, bubble pressure `Pdot`, velocity field, fluxes) is carried as a pair of its value and its sparse derivative with respect to the state, and the chain rule is applied term by term.
- The bubble-wall temperature enters through the wall balance: the exact quadratic root `f_wall_theta_exact` without mass transfer, whose derivative covers every node in the wall rows (three per field for finite differences, all of them for Chebyshev), and the implicit function theorem on `f_wall_mass_bc` (`d(theta_w)/dx = -(dr/dx)/(dr/d(theta_w))`) with mass transfer.
- The low-dimensional closure (radial equation, stress model, non-Newtonian viscosity, medium dissipation) is differentiated by central differences of the closure alone and chained through `Pdot`.

## Inputs
- `t`, `X`: time and state `[R; Rdot; P; theta; Tm; kv; Z]`.
- `jp`: struct set up by `f_imr_nodal` with the state indices, differentiation matrices, grids, dimensionless constants, the wall-balance constants, and handles to the closure, the medium dissipation and the wall solve.

## Verification
Against a central-difference Jacobian of the right-hand side over radial 1–4, stress 0/1/3/5, viscosity models 0/1/5, with and without medium heat and mass transfer, the maximum scaled difference is below `3e-9` for finite differences and `4e-10` for Chebyshev (mass transfer is finite difference only).
