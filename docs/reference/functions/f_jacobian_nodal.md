# `f_jacobian_nodal`

`f_jacobian_nodal(t, X, jp)` returns the analytic Jacobian of the nodal right-hand side, for either differentiation scheme, for `bubtherm = 1`, with or without medium heat transfer (`medtherm`) and vapor mass transfer (`masstrans`). `f_imr_nodal` passes it to the stiff integrators when `jacobian = 2`.

## Inputs

- `t`, `X`: time and state `[R; Rdot; P; theta; Tm; kv; Z]`.
- `jp`: struct set up by `f_imr_nodal` with the state indices, differentiation matrices, grids, dimensionless constants, the wall-balance constants, and handles to the closure, the medium dissipation and the wall solve.

For the algorithm and assumptions, see the [topic guide](../../method/time-integration.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
