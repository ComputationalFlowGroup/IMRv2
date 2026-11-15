# `imr.approximate_initial_stress`

`imr.approximate_initial_stress` computes an approximate initial Cauchy stress at the bubble maximum radius. It uses analytical corrections for viscoelastic, compressibility, and compressive effects.
This estimate is faster than the numerical shooting method, but it is less accurate. The solvers do not call this function. The repository keeps it as a standalone estimate.

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

- This function throws an error if the aggregated stress correction (`fsum`) exceeds 1.
- This analytical method is faster than `imr.physics.initial_stress`'s numerical shooting method, but it is less accurate.

For the algorithm and assumptions, see the [topic guide](../../method/initial-stress.md). See [historical reports](../../verification/convergence.md) for the numerical measurements.

This standalone helper calls `hypergeom`. Check that your MATLAB license includes the needed toolbox. Neither solver calls this helper.
