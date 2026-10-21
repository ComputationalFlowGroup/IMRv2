# `f_viscosity`

`f_viscosity` computes the non-Newtonian viscosity integrals that enter the stress integral of the radial equations. The viscosity is `mu = mu8 + (mu0 - mu8)*f(x)` with `x = lambda*|gammadot|`, and each model supplies `f(x)`.

## Outputs

| Name     | Description |
|----------|-------------|
| `f`      | Viscosity factor at the bubble wall |
| `intf`   | `-I1`; enters the stress integral as `-6*intf*iDRe` (see [`f_stress`](f_stress.md)) |
| `dintf`  | `Rdot` part of `d(intf)/dt`, `-d(gammadot_R*R^3)/dt*I2 + h` with the Leibniz term `h = f*gammadot_R*Rdot/R` |
| `ddintf` | `Rddot` coefficient of `d(intf)/dt`, `2*R^2*I2` |
| `fq`     | Optional: viscosity factor `f(x_R*yq^-3)` at the radii `r/R = yq` (seventh input). The solvers pass the medium grid and use it for the viscous dissipation. |

## Notes

- `R` and `Rdot` may be row vectors (the batched right-hand side of [`f_imr_nodal`](../../method/time-integration.md#batched-right-hand-side)), giving one output column per state. With one state the result is the evaluated in the same elementwise order.

For the algorithm and assumptions, see the [topic guide](../../theory/viscosity-integrals.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
