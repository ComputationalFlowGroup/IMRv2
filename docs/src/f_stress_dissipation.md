# `f_stress_dissipation.m` Documentation

## Overview
`f_stress_dissipation` computes the mechanical dissipation `τ:∇u` that heats the medium when `medtherm = 1`. Both solvers call it.

## Major Sections
1. **Stress-field dissipation**: with the spectral stress representation (`spectral = 1`, `f_imr_spectral` only) the dissipation comes from the stress field itself.
2. **No stress**: zero for stress 0.
3. **Closed-form dissipation**: otherwise the radial stretch terms `Rst`, `x2`, `ix2` and `x4` are computed, followed by the neo-Hookean form (stress 1, 3 and 5) or the quadratic neo-Hookean form (stress 2 and 4). For stress 3–5 this is the Kelvin–Voigt form, which is inconsistent for history-dependent stress; `f_checkers` warns (`IMR:dissipation`).

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
| `fnu`, `DRe`   | Non-Newtonian viscosity factor on the medium grid (from the `fq` output of `f_viscosity`; scalar 0 for Newtonian runs), and the non-Newtonian Reynolds number | –     |

## Outputs

| Name        | Description                           | Units |
|-------------|---------------------------------------|-------|
| `taugradu`  | Mechanical dissipation (`τ:∇u`)        | –     |

## Notes
- The viscous part uses the local viscosity `visc = 1/Re8 + fnu/DRe` (`1/Re8` when `DRe = 0`), with `fnu` evaluated at the local shear rate `gammadot_R*(R/r)^3` of each medium grid point.
- The stress-field branch uses `ZZT*(X(ivisco1)-X(ivisco2))`, the difference of the `tau_rr` and `tau_tt` fields on the medium grid.
- `R`, `Rdot` and the other state inputs may be row vectors, with `X` holding one state per column (the batched right-hand side of [`f_imr_nodal`](f_imr_nodal.md#batched-right-hand-side)). Every operation is elementwise across the columns, in the same order as for one state, so a single state gives the same result as before bit for bit. The output then has one column per state.
