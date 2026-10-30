# `imr.physics.material_stress`

`imr.physics.material_stress` computes the stress integral `S` and its time derivative for the closed-form and single-state viscoelastic models (Kelvin–Voigt, linear Kelvin–Voigt, Maxwell, Jeffreys, Zener, UCM, Oldroyd-B). Both solvers call it on every right-hand-side evaluation.

## Inputs

| Name         | Description                            | Units |
|--------------|----------------------------------------|-------|
| `stress`     | Model selector (0–5, 8); stress 6, 7 (PTT, Giesekus) never reach this function, since `imr.configuration.validate_case` forces the spectral stress field for them | –     |
| `X`          | State vector, which holds the stress states | –     |
| `Req`, `R`   | Equilibrium and current radius         | –     |
| `Ca`, `De`   | Cauchy and Deborah numbers             | –     |
| `Re8`        | Reynolds number                        | –     |
| `Rdot`       | Wall velocity                          | –     |
| `alphax`     | Strain-stiffening parameter of the quadratic models | –     |
| `ivisco1`    | Index of the first stress state in `X`  | –     |
| `ivisco2`    | Index of the second stress state in `X` | –     |
| `LAM`        | Ratio of retardation to relaxation time, `lambda2/lambda1` | –     |
| `intfnu`, `dintfnu` | Non-Newtonian stress integral and the `Rdot` part of its time derivative, from `imr.physics.viscosity` | –     |
| `iDRe`       | Scaling of the non-Newtonian stress: `1/DRe` for the memoryless Kelvin-Voigt family (stress 1, 2, 8) and `LAM/DRe` for the Maxwell/Zener and UCM family (stress 3–5) (0 when `DRe = 0`), set by `imr.configuration.prepare_case` | –     |

## Outputs

| Name    | Description                        | Units |
|---------|------------------------------------|-------|
| `S`     | Stress integral                    | –     |
| `Sdot`  | Time derivative of `S`             | –     |
| `Z1dot` | Time derivative of the first stress state (stress 3–5), empty otherwise  | –     |
| `Z2dot` | Time derivative of the second stress state (stress 5), empty otherwise | –     |

## Notes

- The spectral stress representation (`imr.solve_spectral`, `stress = 5, 6, 7`, `spectral = 1`) does not use this function; `f_spectral_stress` in `imr.solvers.build_spectral` implements it instead, including the PTT and Giesekus closures (see [stress and viscosity](../../model/stress-viscosity.md#ptt-and-giesekus-stress-6-7)).
- The non-Newtonian contribution enters every stress model with a viscous term (stress 1–5, 8) as `-6*intfnu*iDRe` in `S` and `-6*dintfnu*iDRe` in `Sdot`, so that `Sdot` is the time derivative of `S` (the `Rddot` part is carried by `ddintfnu` in [`imr.physics.radial_acceleration`](f_radial_eq.md)). For stress 3–5 the non-Newtonian viscosity takes the place of the solvent viscosity `LAM/Re8`, which gives the `LAM/DRe` scaling; stress 8 uses `1/DRe` directly, like stress 1, 2. `imr.configuration.validate_case` rejects `nu_model > 0` with stress 0.
- Stress 8 (linear Kelvin-Voigt, Yang-Church) has the same closed-form, no-history structure as stress 1, 2, with a spring linear in `Rst = Req/R` instead of a finite-strain invariant; see [stress and viscosity](../../model/stress-viscosity.md#linear-kelvin-voigt-stress-8).
- `R`, `Rdot` and the other state inputs may be row vectors, with `X` holding one state per column (the batched right-hand side of [`imr.solve_nodal`](../../method/time-integration.md#batched-right-hand-side)). Every operation is elementwise across the columns, in the same order as for one state, so a single state gives the same result as before bit for bit.
