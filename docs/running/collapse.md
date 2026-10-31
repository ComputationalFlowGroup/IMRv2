# Collapse initial conditions

`collapse = 1` initializes the bubble at maximum radius with zero wall velocity. `imr.configuration.prepare_case` constructs that state before time integration; it is not a separate radial equation.

Use [compatibility](../reference/compatibility.md) for the exact prerequisites. The current checker requires bubble heat, medium heat, vapor pressure and mass transfer, so supported collapse uses the nodal solver with radial 1 or 2.

Preparation solves the initial gas/vapor/surface-tension balance for a positive equilibrium radius $r$:

$$
(1-P_v^*)r^3+\frac{r^2}{\mathrm{We}}-P_{g0}^*=0.
$$

Here $P_{g0}^*$ is `p0star`, the gas-pressure coefficient at the initial radius $R=1$, and $P_v^*$ is the enabled saturation pressure scaled by ambient pressure. Dividing by $r^3$ recovers gas plus vapor pressure balanced by ambient pressure and surface tension. Preparation brackets the stable positive root; if no such equilibrium exists, it raises an input error. This replaces an unconstrained solve in vapor mass. The [settings suite](../validation/settings.md) checks the radius against an independently evaluated pressure balance, including zero surface tension and boiling-pressure limits.

| Stress | Initial history | User action |
|---|---|---|
| 0–4 | No stress history state (Kelvin-Voigt family) | Start from rest at maximum radius |
| 5–6 | `imr.physics.initial_stress` shoots a subsonic growth trajectory unless `stress0` is supplied | Provide `stress0` if no admissible growth reaches maximum radius |
| 7–11 | No automatic maximum-radius stress calculation | Supply `stress0` (the two-state form for stress 9 with `spectral = 0`, or `2*Nv` modes otherwise); medium heat transfer remains incompatible in the nodal solver, and stress 10, 11 are rejected there outright |

The stress-9 requirement is necessary but does not remove its medium-heating incompatibility. A full supported collapse configuration with stress 9–11 is therefore not currently available through the checker; do not infer support from the initial-stress rule alone. `imr.physics.initial_stress` is likewise not extended to the 2nd-order Fung and Mooney-Rivlin Zener springs (7, 8), unlike their Kelvin-Voigt variants (3, 4).

For stress 6, the initial stress calculation includes `alphax`; stress 5 uses zero for that correction. The approximate standalone helper is not called by the solver. See [initial-stress method](../method/initial-stress.md) and [troubleshooting](troubleshooting.md).

The recomputed equilibrium radius and zero initial wall velocity are setup quantities. The existing historical report that the recomputed radius reproduces the input over selected cases is preserved in [reported comparisons](../validation/convergence.md#reported-collapse-setup-comparison), not presented as a new run.
