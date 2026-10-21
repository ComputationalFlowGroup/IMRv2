# `f_imr_spectral`

`f_imr_spectral` integrates bubble dynamics and thermal fields with Chebyshev expansions. It also provides a spectral representation of UCM/Oldroyd-B polymer stresses.

```matlab
[t, R, Rdot, P, T, Tm, kv, inputs] = f_imr_spectral(varargin);
```

Inputs follow the same dimensionless [configuration forms](../running/configuration.md) as the nodal solver. The [parameter table](../model/parameters.md) owns accepted names and defaults. `kv`, output 7, is always empty; other arrays follow the [shared output contract](../running/output.md).

| Input | Spectral-specific behavior |
|---|---|
| `spectral` | Selects stress-5 representation, not the solver function |
| `nv`, `lv` | Stress expansion size and mapping length |
| `nt`, `mt`, `lt` | Thermal expansion sizes and medium mapping length |
| `jacobian = 2` | Default hybrid analytic Jacobian |
| `jacobian = 0` or `1` | Dense numerical Jacobian |
| `jacobian = 3` | Rejected; nodal solver only |
| `scheme` | Leave at zero; nonzero value is rejected |

Vapor mass transfer and currently supported collapse are unavailable. With stress 5 and medium heat, the checker selects the local stress field and warns. See [compatibility](../reference/compatibility.md).

Implementation: `src/f_imr_spectral.m`, `f_call_params`, `f_odesolve`. See [discretization](../method/discretization.md), [stress heating theory](../theory/stress-heating.md) and [reported comparisons](../validation/convergence.md).
