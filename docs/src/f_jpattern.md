# `f_jpattern.m` Documentation

## Overview
`Sp = f_jpattern(n, ith, itm, ikv, iz)` builds the Jacobian sparsity pattern of the finite difference solver from its state indices. `f_imr_fd` passes it to the stiff integrators when `jacobian = 1`, and when `jacobian = 2` without bubble heat transfer, where no analytic Jacobian exists.

## Inputs
- `n`: length of the state vector `[R, Rdot, P, theta, Tm, kv, Z]`.
- `ith`, `itm`, `ikv`, `iz`: indices of the bubble temperature, medium temperature, vapor concentration and stress states.

## Notes
- The pattern is a conservative superset of the true Jacobian. The bubble pressure equation and the wall boundary condition feed every field row, so all field rows couple to `R`, `Rdot`, `P` and the wall stencil columns, while the field operators themselves are tridiagonal.
- GNU Octave's `ode15s` ignores `JPattern`.
