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
| `iDRe`       | `1/DRe`, the inverse non-Newtonian Reynolds number | –     |

## Outputs

| Name    | Description                        | Units |
|---------|------------------------------------|-------|
| `S`     | Stress integral                    | –     |
| `Sdot`  | Time derivative of `S`             | –     |
| `Z1dot` | Time derivative of the first stress state (stress 3–5), empty otherwise  | –     |
| `Z2dot` | Time derivative of the second stress state (stress 5), empty otherwise | –     |

## Notes
- The spectral stress representation (`f_imr_spectral`, `stress = 5`, `spectral = 1`) does not use this function.
- The non-Newtonian contribution enters as `-6*intfnu*iDRe`.
