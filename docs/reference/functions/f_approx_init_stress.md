# `imr.approximate_initial_stress`

`imr.approximate_initial_stress` computes an approximate initial Cauchy stress at the bubble maximum radius using analytical corrections for viscoelastic, compressibility, and compressive effects.
This provides a lower-fidelity but faster estimate compared to the numerical shooting method. The solvers do not call it; it is kept as a standalone estimate.

## Inputs

| Name     | Description                                       | Units   |
|----------|---------------------------------------------------|---------|
| `Ro`     | Non-dimensional bubble radius (`R/R_eq`)          | –       |
| `kappa`  | Ratio of specific heats of gas                    | –       |
| `al_nd`  | Parameter used as `Mc = 1/al_nd`; physical interpretation needs confirmation         | –       |
| `pwv_nd` | Non-dimensional vapor pressure ratio              | –       |
| `We`     | Weber number (surface tension effects)            | –       |
| `Re`     | Reynolds number (viscous effects)                 | –       |
| `De`     | Deborah number (viscoelastic relaxation effects)  | –       |
| `Ca`     | Cauchy number (elastic to inertial ratio)         | –       |
| `alpha`  | Quadratic KV model exponent                       | –       |

## Outputs

| Name        | Description                             | Units   |
|-------------|-----------------------------------------|---------|
| `Smaxpred`  | Predicted maximum Cauchy stress ratio   | –       |

## Notes

- Throws an error if the aggregated stress correction (`fsum`) exceeds 1.
- Analytical method is faster but less accurate than `imr.physics.initial_stress` numerical shooting.

For the algorithm and assumptions, see the [topic guide](../../method/initial-stress.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).

This standalone helper calls `hypergeom`; verify its toolbox availability. Neither solver calls it.
