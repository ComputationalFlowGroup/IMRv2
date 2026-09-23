# `f_stress.m` Documentation

## Overview
`f_stress` computes the stress integral and its time derivative for a range of viscoelastic models (Kelvin–Voigt, Maxwell, Jeffreys, Zener, Oldroyd-B).

## Major Sections
1. **No Stress (Model 0)**
2. **Kelvin–Voigt Neo-Hookean (Model 1)**
3. **Quadratic KV Neo-Hookean (Model 2)**
4. **Linear Maxwell/Jeffreys/Zener (Model 3)**
5. **Quadratic Maxwell Variant (Model 4)**
6. **Upper-Convected Maxwell / Oldroyd-B (Model 5)**
7. **Special Placeholder (Model -1)**

## Inputs

| Name         | Description                            | Units |
|--------------|----------------------------------------|-------|
| `stress`     | Model selector (0–5)                  | –     |
| `X`          | State vector containing auxiliary vars | –     |
| `Req`, `R`   | Equilibrium and current radius         | –     |
| `Ca`, `De`   | Cauchy and Deborah numbers             | –     |
| `Re8`        | Reynolds number                        | –     |
| `Rdot`       | Wall velocity                          | –     |
| `alphax`     | Quadratic KV exponent                  | –     |
| `ivisco1`    | Index of first viscous state variable  | –     |
| `ivisco2`    | Index of second viscous state variable | –     |
| `LAM`        | Dimensionless elastic modulus          | –     |
| `intfnu`, `dintfnu` | Non-Newtonian integrals        | –     |
| `iDRe`       | Inverse Debye–Stokes number           | –     |

## Outputs

| Name    | Description                        | Units |
|---------|------------------------------------|-------|
| `S`     | Stress integral                    | –     |
| `Sdot`  | Time derivative of `S`             | –     |
| `Z1dot` | First auxiliary stress derivative  | –     |
| `Z2dot` | Second auxiliary stress derivative | –     |

## Notes
- The spectral stress representation (`f_imr_spectral`, `stress = 5`, `spectral = 1`) does not use this function.
- The placeholder model (-1) is for internal code consistency.
