# `imr.physics.viscosity`

`imr.physics.viscosity` computes the non-Newtonian viscosity integrals that enter the stress integral of the radial equations. The viscosity is $\mu=\mu_8+(\mu_0-\mu_8)f(x)$ with $x=\lambda|\dot\gamma|$, and each model supplies $f(x)$.

## Outputs

| Name     | Description |
|----------|-------------|
| `f`      | Viscosity factor at the bubble wall |
| `intf`   | $-I_1$; enters the stress integral as $-6\,\mathtt{intf}\,\mathtt{iDRe}$ (see [`imr.physics.material_stress`](material_stress.md)) |
| `dintf`  | $\dot R$ part of $d(\mathtt{intf})/dt$, $-\dfrac{d(\dot\gamma_R R^3)}{dt}I_2+h$ with the Leibniz term $h=f\dot\gamma_R\dot R/R$ |
| `ddintf` | $\ddot R$ coefficient of $d(\mathtt{intf})/dt$, $2R^2I_2$ |
| `fq`     | Optional: viscosity factor $f(x_R y_q^{-3})$ at the radii $r/R=y_q$ (seventh input). The solvers pass the medium grid and use it for the viscous dissipation. |

## Notes

- `R` and `Rdot` can be row vectors (the batched right-hand side of [`imr.solve_nodal`](../../method/time-integration.md#batched-right-hand-side)). This gives one output column per state. With one state, the function evaluates the result in the same elementwise order.

For the algorithm and assumptions, see the [topic guide](../../model/stress-viscosity.md#non-newtonian-viscosity). Numerical measurements are retained as [historical reports](../../verification/convergence.md).
