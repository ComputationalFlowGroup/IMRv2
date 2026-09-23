# `f_jacobian_spectral.m` Documentation

## Overview
`f_jacobian_spectral(t, X, jp)` returns the analytic Jacobian of the spectral right-hand side, passed by `f_imr_spectral` to the stiff integrators when `jacobian = 2`.

## Method
- Quantities are carried as value/derivative pairs (`f_dual_ops`), so the chain rule through the bubble temperature, the bubble pressure, the medium temperature and the spectral stress fields is exact.
- The thermal boundary condition solve `A(s)*q = rhs` with `A(s) = A0 - (s-1)*e1*v'` (`s` depends on the wall temperature) is differentiated as `dq = A(s)\(d(rhs) + e1*(v'*q)*ds)`, reusing the prefactored `A0` and the Sherman–Morrison update.
- The radial equation and the closed-form stress models are differentiated by central differences of the low-dimensional closure alone.

## Verification
Against a central-difference Jacobian over radial 1–4, stress 0/1/3/5, Newtonian and Carreau viscosity, polytropic, cold-medium and full thermal models, and the spectral stress representation (102 configurations), the maximum scaled difference is `3e-7`.
