# Collapse initial conditions

`collapse = 1` initializes the bubble at maximum radius, with zero wall velocity. `imr.configuration.prepare_case` constructs that state before time integration. This is not a separate radial equation.

Use [compatibility](../reference/compatibility.md) for the exact requirements. The current checker requires bubble heat, medium heat, vapor pressure and mass transfer. So supported collapse uses the nodal solver with radial 1 or 2.

Preparation solves the initial gas/vapor/surface-tension balance for a positive equilibrium radius $r$:

$$
(1-P_v^*)r^3+\frac{r^2}{\mathrm{We}}-P_{g0}^*=0.
$$

Here $P_{g0}^*$ is `p0star`, the gas-pressure coefficient at the initial radius $R=1$. $P_v^*$ is the enabled saturation pressure, scaled by ambient pressure. Dividing by $r^3$ recovers gas plus vapor pressure, balanced by ambient pressure and surface tension. Preparation brackets the stable positive root. If no such equilibrium exists, preparation raises an input error. This replaces an unconstrained solve in vapor mass. The [settings suite](../verification/settings.md) checks the radius against an independently evaluated pressure balance. This check includes the zero surface tension limit and the boiling-pressure limit.

| Stress | Initial history | User action |
|---|---|---|
| 0–4 | No stress history state (Kelvin-Voigt family) | Start from rest at maximum radius |
| 5–6 | `imr.physics.initial_stress` shoots a subsonic growth trajectory unless `stress0` is supplied | Provide `stress0` if no admissible growth reaches maximum radius |
| 7–11 | No automatic maximum-radius stress calculation | Supply `stress0` (the two-state form for stress 9 with `spectral = 0`, or `2*Nv` modes otherwise). Medium heat transfer stays incompatible in the nodal solver. The nodal solver rejects stress 10 and 11 outright |

The stress-9 requirement is necessary, but it does not remove the medium-heating incompatibility. So the checker does not currently support a full collapse configuration with stress 9–11. Do not infer support from the initial-stress rule alone. `imr.physics.initial_stress` is also not extended to the 2nd-order Fung and Mooney-Rivlin Zener springs (7, 8). Their Kelvin-Voigt variants (3, 4) do have this extension.

For stress 6, the initial stress calculation includes `alphax`. For stress 5, the calculation uses zero for that correction. The solver does not call the approximate standalone helper. See [initial-stress method](../method/initial-stress.md) and [troubleshooting](troubleshooting.md).

The recomputed equilibrium radius and the zero initial wall velocity are setup quantities. [Reported comparisons](../verification/convergence.md#reported-collapse-setup-comparison) keeps an existing historical report. That report states that the recomputed radius reproduces the input over selected cases. The report is historical. It is not a new run.
