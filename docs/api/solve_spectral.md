# `imr.solve_spectral`

`imr.solve_spectral` integrates bubble dynamics and thermal fields with Chebyshev expansions. It also provides a spectral representation of UCM/Oldroyd-B, PTT and Giesekus polymer stresses.

```matlab
[t, R, Rdot, P, T, Tm, kv, inputs, Rddot, S, Sdot] = imr.solve_spectral(varargin);
```

Inputs follow the same dimensionless [configuration forms](../running/configuration.md) as the nodal solver. The [parameter table](../model/parameters.md) owns accepted names and defaults. `kv`, output 7, is always empty; other arrays follow the [shared output contract](../running/output.md).

The first eight outputs retain their positions; optional output 9 is wall acceleration, and optional outputs 10 and 11 are the stress integral and its time derivative. Use [`imr.solve(...,'solver','spectral')`](solve.md) for the common entry point. This wrapper fixes the spectral solver and rejects a conflicting `solver` value.

| Input | Spectral-specific behavior |
|---|---|
| `spectral` | Selects the stress-9 field representation, not the solver function; always on for stress 10 (PTT) and 11 (Giesekus), which have no other representation |
| `nv`, `lv` | Stress expansion size and mapping length |
| `nt`, `mt`, `lt` | Thermal expansion sizes and medium mapping length |
| `jacobian = 2` | Default hybrid analytic Jacobian |
| `jacobian = 0` or `1` | Dense numerical Jacobian |
| `jacobian = 3` | Rejected; nodal solver only |
| `scheme` | Leave at zero; nonzero value is rejected |

Vapor mass transfer and currently supported collapse are unavailable. With stress 9 and medium heat, the checker selects the local stress field and warns. Stress 10 and 11 (PTT, Giesekus; see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-10-11)) always use it, forcing `spectral = 1` with a warning if it was 0. See [compatibility](../reference/compatibility.md).

Implementation: `imr.configuration.prepare_case` → `imr.solvers.build` → `imr.numerics.integrate` → output reconstruction. Spectral operators and closures live in `src/+imr/+solvers/build_spectral.m`. See [discretization](../method/discretization.md), [stress heating theory](../theory/stress-heating.md) and [reported comparisons](../validation/convergence.md).
