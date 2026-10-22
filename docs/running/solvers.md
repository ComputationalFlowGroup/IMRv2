# Choose a solver

IMR has two solver functions and three spatial modes. Select a mode for the fields you need, then verify resolution for your case.

| Mode | Invocation | Main distinction |
|---|---|---|
| Nodal finite differences | `imr.solve_nodal(...,'scheme',0)` | Default, second-order spatial differences |
| Nodal Chebyshev | `imr.solve_nodal(...,'scheme',1)` | Same physical code with Chebyshev differentiation matrices |
| Spectral | `imr.solve_spectral(...)` | Chebyshev expansions, including optional UCM/Oldroyd-B stress field |

Both functions use [`imr.numerics.integrate`](../reference/functions/f_odesolve.md) for time integration. `scheme` selects only the nodal differentiation matrices; the input named `spectral` selects a stress representation, not the solver function.

Use the [compatibility matrix](../reference/compatibility.md) for physical restrictions. In particular, vapor mass transfer and supported collapse initialization require the nodal solver. Stress 5 with medium heat transfer needs the spectral solver's local stress field. Radial 5–6 support only a polytropic bubble.

See [discretization](../method/discretization.md), [integration and Jacobians](../method/time-integration.md), and [reported convergence](../validation/convergence.md) before treating a grid as resolved. Neither the default grid nor an integrator tolerance alone establishes spatial accuracy.
