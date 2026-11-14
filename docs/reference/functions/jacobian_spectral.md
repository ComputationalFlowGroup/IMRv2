# `imr.numerics.jacobian_spectral`

`imr.numerics.jacobian_spectral(t, X, jp)` returns the analytic Jacobian of the spectral right-hand side. `imr.solve_spectral` passes this Jacobian to the stiff integrators when `jacobian = 2`.

For the algorithm and assumptions, see the [topic guide](../../method/time-integration.md). See [historical reports](../../verification/convergence.md) for the numerical measurements.
