# `f_stress_dissipation.m` Documentation

## Overview
`f_stress_dissipation` computes the mechanical dissipation `τ:∇u` that heats the medium when `medtherm = 1`. Both solvers call it.

## Major Sections
1. **Radial Stretch Calculation**: Compute `Rst`, `x2`, `ix2`, `x4`.
2. **Closed-form dissipation**: one branch per stress model 0–5. For stress 3–5 this is the Kelvin–Voigt form, which is inconsistent for history-dependent stress; `f_checkers` warns (`IMR:dissipation`).
3. **Stress-field dissipation**: with the spectral stress representation (`spectral = 1`, `f_imr_spectral` only) the dissipation comes from the stress field itself and replaces the closed form.

## Inputs

| Name           | Description                            | Units |
|----------------|----------------------------------------|-------|
| `stress`       | Stress model selector (0–5)            | –     |
| `spectral`     | Spectral stress representation flag (1 only in `f_imr_spectral` with stress 5) | –     |
| `Req`, `R`     | Equilibrium and current radius         | –     |
| `Rdot`         | Wall velocity                          | –     |
| `Ca`, `Br`     | Cauchy number and Brinkman number      | –     |
| `Re8`, `alphax`| Reynolds number and strain-stiffening parameter | –     |
| `yT2`, `yT3`   | Powers of the medium radial coordinate `r/R` | –     |
| `iyT3`, `iyT4`, `iyT6` | Inverses of grid powers        | –     |
| `X`, `ZZT`     | State vector and the matrix from stress coefficients to medium grid values | –     |
| `ivisco1`, `ivisco2` | Indices of the stress coefficients in `X` | –     |
| `fnu`, `DRe`   | Non-Newtonian viscosity factor at the wall, and the non-Newtonian Reynolds number | –     |

## Outputs

| Name        | Description                           | Units |
|-------------|---------------------------------------|-------|
| `taugradu`  | Mechanical dissipation (`τ:∇u`)        | –     |

## Notes
- The stress-field branch uses `ZZT*(X(ivisco1)-X(ivisco2))`, the difference of the `tau_rr` and `tau_tt` fields on the medium grid.
