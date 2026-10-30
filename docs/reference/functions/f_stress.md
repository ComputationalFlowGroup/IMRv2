# `imr.physics.material_stress`

`imr.physics.material_stress` computes the stress integral `S` and its time derivative for the closed-form and single-state viscoelastic models (Kelvin–Voigt, Maxwell, Jeffreys, Zener, UCM, Oldroyd-B, and the 2nd-order Fung and Mooney-Rivlin Kelvin-Voigt/Zener springs). Both solvers call it on every right-hand-side evaluation.

## Inputs

| Name         | Description                            | Units |
|--------------|----------------------------------------|-------|
| `stress`     | Model selector (0–5, 8–11); stress 6, 7 (PTT, Giesekus) never reach this function, since `imr.configuration.validate_case` forces the spectral stress field for them | –     |
| `X`          | State vector, which holds the stress states | –     |
| `Req`, `R`   | Equilibrium and current radius         | –     |
| `Ca`, `De`   | Cauchy and Deborah numbers             | –     |
| `Re8`        | Reynolds number                        | –     |
| `Rdot`       | Wall velocity                          | –     |
| `alphax`     | Strain-stiffening parameter of the quadratic models (2, 4); the next series term reuses it for the 2nd-order Fung spring (8, 9) | –     |
| `ivisco1`    | Index of the first stress state in `X`  | –     |
| `ivisco2`    | Index of the second stress state in `X` | –     |
| `LAM`        | Ratio of retardation to relaxation time, `lambda2/lambda1` | –     |
| `intfnu`, `dintfnu` | Non-Newtonian stress integral and the `Rdot` part of its time derivative, from `imr.physics.viscosity` | –     |
| `iDRe`       | Scaling of the non-Newtonian stress: `1/DRe` for the Kelvin-Voigt family (1, 2, 8, 10) and `LAM/DRe` for stress 3–5, 9, 11 (0 when `DRe = 0`), set by `imr.configuration.prepare_case` | –     |
| `C2`         | Mooney-Rivlin second-modulus ratio (10, 11); unused otherwise | –     |

## Outputs

| Name    | Description                        | Units |
|---------|------------------------------------|-------|
| `S`     | Stress integral                    | –     |
| `Sdot`  | Time derivative of `S`             | –     |
| `Z1dot` | Time derivative of the first stress state (stress 3–5, 9, 11), empty otherwise  | –     |
| `Z2dot` | Time derivative of the second stress state (stress 5), empty otherwise | –     |

## Notes

- The spectral stress representation (`imr.solve_spectral`, `stress = 5, 6, 7`, `spectral = 1`) does not use this function; `f_spectral_stress` in `imr.solvers.build_spectral` implements it instead, including the PTT and Giesekus closures (see [stress and viscosity](../../model/stress-viscosity.md#ptt-and-giesekus-stress-6-7)).
- The non-Newtonian contribution enters every stress model with a viscous term (stress 1–5, 8–11) as `-6*intfnu*iDRe` in `S` and `-6*dintfnu*iDRe` in `Sdot`, so that `Sdot` is the time derivative of `S` (the `Rddot` part is carried by `ddintfnu` in [`imr.physics.radial_acceleration`](f_radial_eq.md)). For stress 3–5, 9, 11 the non-Newtonian viscosity takes the place of the solvent viscosity `LAM/Re8`, which gives the `LAM/DRe` scaling. `imr.configuration.validate_case` rejects `nu_model > 0` with stress 0.
- The Zener variants of the new springs (stress 9, 11) set `Ze = R^3` times the corresponding Kelvin-Voigt spring's own elastic-only `S` (8, 10), the same relation stress 3 and 4's `Ze` already has to stress 1 and 2's elastic `S`.
- `R`, `Rdot` and the other state inputs may be row vectors, with `X` holding one state per column (the batched right-hand side of [`imr.solve_nodal`](../../method/time-integration.md#batched-right-hand-side)). Every operation is elementwise across the columns, in the same order as for one state, so a single state gives the same result as before bit for bit.
