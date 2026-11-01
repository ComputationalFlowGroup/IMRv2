# `imr.numerics.integrate`

`imr.numerics.integrate` is a wrapper that selects and executes a stiff ODE solver (`ode15s`, `ode23tb` or CVODES) with user-defined time stepping and tolerance options.

## Inputs

| Name       | Description                                  | Units   |
|------------|----------------------------------------------|---------|
| `bubble`   | Function handle for ODE system               | –       |
| `init`     | Initial condition vector                     | –       |
| `method`   | Solver choice (3, 15 or 23)                  | –       |
| `divisions`| Minimum number of time steps (0 = automatic) | –       |
| `tspan`    | Time span or output time vector              | –       |
| `Sp`       | Optional: Jacobian sparsity pattern (`JPattern`) or Jacobian function handle | – |
| `rtol`     | Optional: relative tolerance (default `1e-8`) | –      |
| `atol`     | Optional: absolute tolerance (default `1e-8`) | –      |

## Outputs

| Name | Description              | Units   |
|------|--------------------------|---------|
| `t`  | Time vector              | –       |
| `X`  | Solution matrix (states) | –       |

## Notes

- Throws an error for unsupported `method` codes.
- `MaxStep` is set to `(tspan(end)-tspan(1))/divisions` when `divisions > 0`.
- CVODES: the `ode` object's properties are set by assignment (`F.Solver = 'cvodesstiff'` and so on), so the file also parses where the object is absent, and a clear error names the other methods there. A Jacobian function becomes `odeJacobian('Jacobian', Sp)` and a pattern `odeJacobian('SparsityPattern', Sp)`, and `MaxStep` goes to `F.SolverOptions`. With two times in `tspan` the solver's own steps are returned, otherwise the solution at the times of `tspan`, as with `ode15s` and `ode23tb`.
- CVODES is compiled, but it still calls the right-hand side and the Jacobian in MATLAB, so it saves the integrator's own interpreted work (step control, error estimates, the Newton iteration), not the model evaluations. `benchmarks/integrators/s_compare_integrators.m` measures what that is worth.

For the algorithm and assumptions, see the [topic guide](../../method/time-integration.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
