# `f_stress.m` Documentation

## Overview
`f_stress` computes the stress integral `S` and its time derivative for the closed-form and single-state viscoelastic models (Kelvin–Voigt, Maxwell, Jeffreys, Zener, UCM, Oldroyd-B). Both solvers call it on every right-hand-side evaluation.

## Major Sections
1. **No Stress (Model 0)**
2. **Kelvin–Voigt Neo-Hookean (Model 1)**
3. **Quadratic KV Neo-Hookean (Model 2)**
4. **Linear Maxwell/Jeffreys/Zener (Model 3)**
5. **Quadratic Maxwell Variant (Model 4)**
6. **Upper-Convected Maxwell / Oldroyd-B (Model 5)**, through the two stress states `Z1`, `Z2`

## Inputs

| Name         | Description                            | Units |
|--------------|----------------------------------------|-------|
| `stress`     | Model selector (0–5)                  | –     |
| `X`          | State vector, which holds the stress states | –     |
| `Req`, `R`   | Equilibrium and current radius         | –     |
| `Ca`, `De`   | Cauchy and Deborah numbers             | –     |
| `Re8`        | Reynolds number                        | –     |
| `Rdot`       | Wall velocity                          | –     |
| `alphax`     | Strain-stiffening parameter of the quadratic models | –     |
| `ivisco1`    | Index of the first stress state in `X`  | –     |
| `ivisco2`    | Index of the second stress state in `X` | –     |
| `LAM`        | Ratio of retardation to relaxation time, `lambda2/lambda1` | –     |
| `intfnu`, `dintfnu` | Non-Newtonian stress integral and the `Rdot` part of its time derivative, from `f_viscosity` | –     |
| `iDRe`       | Scaling of the non-Newtonian stress: `1/DRe` for stress 1–2 and `LAM/DRe` for stress 3–5 (0 when `DRe = 0`), set by `f_call_params` | –     |

## Outputs

| Name    | Description                        | Units |
|---------|------------------------------------|-------|
| `S`     | Stress integral                    | –     |
| `Sdot`  | Time derivative of `S`             | –     |
| `Z1dot` | Time derivative of the first stress state (stress 3–5), empty otherwise  | –     |
| `Z2dot` | Time derivative of the second stress state (stress 5), empty otherwise | –     |

## Notes
- The spectral stress representation (`f_imr_spectral`, `stress = 5`, `spectral = 1`) does not use this function.
- The non-Newtonian contribution enters every stress model with a viscous term (stress 1–5) as `-6*intfnu*iDRe` in `S` and `-6*dintfnu*iDRe` in `Sdot`, so that `Sdot` is the time derivative of `S` (the `Rddot` part is carried by `ddintfnu` in [`f_radial_eq`](f_radial_eq.md)). For stress 3–5 the non-Newtonian viscosity takes the place of the solvent viscosity `LAM/Re8`, which gives the `LAM/DRe` scaling. `f_checkers` rejects `nu_model > 0` with stress 0.
- `R`, `Rdot` and the other state inputs may be row vectors, with `X` holding one state per column (the batched right-hand side of [`f_imr_nodal`](f_imr_nodal.md#batched-right-hand-side)). Every operation is elementwise across the columns, in the same order as for one state, so a single state gives the same result as before bit for bit.
