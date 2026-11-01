# `imr.numerics.jacobian_spectral`

`imr.numerics.jacobian_spectral(t, X, jp)` returns the analytic Jacobian of the spectral right-hand side, passed by `imr.solve_spectral` to the stiff integrators when `jacobian = 2`.

For the algorithm and assumptions, see the [topic guide](../../method/time-integration.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
