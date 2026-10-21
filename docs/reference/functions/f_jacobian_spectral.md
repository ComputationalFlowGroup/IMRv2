# `f_jacobian_spectral`

`f_jacobian_spectral(t, X, jp)` returns the analytic Jacobian of the spectral right-hand side, passed by `f_imr_spectral` to the stiff integrators when `jacobian = 2`.

For the algorithm and assumptions, see the [topic guide](../../method/time-integration.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
